//! The rasterizer on the graphics card, by way of wgpu: Metal on macOS,
//! Vulkan or OpenGL on Linux, DirectX 12 or Vulkan on Windows.
//!
//! GPUI takes images, not draw calls, so the frame is drawn offscreen and read
//! back. It is the same picture as the CPU path draws, pixel for pixel near
//! enough: the shaders in `gpu.wgsl` follow `render.rs` step by step. Where no
//! usable device is found the CPU path draws instead.

use std::sync::{Mutex, OnceLock};

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::render::{Camera, Frame, Scene, Style, OPAQUE};

/// A sphere or a capsule, as the impostor shader takes it. A sphere is a
/// capsule whose ends coincide, told apart so it can keep its ink rim.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Impostor {
    a: [f32; 3],
    r: f32,
    b: [f32; 3],
    kind: f32,
    /// Colour at `a`, and the opacity.
    ca: [f32; 4],
    cb: [f32; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Vertex {
    p: [f32; 3],
    n: [f32; 3],
    c: [f32; 3],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Uniforms {
    rot: [[f32; 4]; 3],
    center: [f32; 4],
    view: [f32; 4],
    size: [f32; 4],
    depth: [f32; 4],
    page: [f32; 4],
    params: [f32; 4],
}

const COLOR: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;
/// Read back by the ink pass, so it doubles as the frame's depth map.
const DEPTH: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
/// What GPUI uploads, so the readback needs no swizzle.
const OUT: wgpu::TextureFormat = wgpu::TextureFormat::Bgra8Unorm;
/// Where the ink lines fall, before they are thickened.
const EDGE: wgpu::TextureFormat = wgpu::TextureFormat::R8Unorm;

pub struct Gpu {
    device: wgpu::Device,
    queue: wgpu::Queue,
    max_size: u32,
    uniform_layout: wgpu::BindGroupLayout,
    post_layout: wgpu::BindGroupLayout,
    edge_layout: wgpu::BindGroupLayout,
    impostors: wgpu::RenderPipeline,
    glass: wgpu::RenderPipeline,
    ribbons: wgpu::RenderPipeline,
    edges: wgpu::RenderPipeline,
    post: wgpu::RenderPipeline,
    uniforms: wgpu::Buffer,
    state: Mutex<State>,
}

/// What is kept between frames: the scene's buffers, rebuilt when the scene
/// changes, and the targets, rebuilt when the size does.
#[derive(Default)]
struct State {
    scene: u64,
    opaque: Option<(wgpu::Buffer, u32)>,
    ribbons: Option<(wgpu::Buffer, u32)>,
    /// Translucent primitives, re-sorted and re-uploaded each frame.
    glass: Vec<Impostor>,
    glass_buffer: Option<(wgpu::Buffer, usize)>,
    targets: Option<Targets>,
}

struct Targets {
    width: u32,
    height: u32,
    color: wgpu::TextureView,
    depth: wgpu::TextureView,
    edge: wgpu::TextureView,
    out: wgpu::Texture,
    out_view: wgpu::TextureView,
    readback: wgpu::Buffer,
    row: u32,
    post: wgpu::BindGroup,
    edges: wgpu::BindGroup,
    uniform: wgpu::BindGroup,
}

/// The device, set up on first use. `None` when there is none to be had, or
/// when `PROFOLD_CPU` is set.
pub fn get() -> Option<&'static Gpu> {
    static GPU: OnceLock<Option<Gpu>> = OnceLock::new();
    GPU.get_or_init(|| {
        if std::env::var_os("PROFOLD_CPU").is_some() {
            return None;
        }
        match pollster::block_on(Gpu::new()) {
            Ok(gpu) => Some(gpu),
            Err(e) => {
                eprintln!("profold: drawing on the CPU ({e})");
                None
            }
        }
    })
    .as_ref()
}

impl Gpu {
    async fn new() -> anyhow::Result<Self> {
        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await?;
        // A software rasterizer such as llvmpipe is no faster than our own.
        if adapter.get_info().device_type == wgpu::DeviceType::Cpu {
            anyhow::bail!("only a software adapter, {}", adapter.get_info().name);
        }
        let limits = wgpu::Limits::downlevel_defaults().using_resolution(adapter.limits());
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("profold"),
                required_limits: limits.clone(),
                ..Default::default()
            })
            .await?;

        let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
        let gpu = Self::build(device, queue, limits.max_texture_dimension_2d);
        if let Some(e) = scope.pop().await {
            anyhow::bail!("{e}");
        }
        Ok(gpu)
    }

    fn build(device: wgpu::Device, queue: wgpu::Queue, max_size: u32) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("profold"),
            source: wgpu::ShaderSource::Wgsl(include_str!("gpu.wgsl").into()),
        });

        let uniform_entry = wgpu::BindGroupLayoutEntry {
            binding: 0,
            visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let texture_entry = |binding, sample_type| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::FRAGMENT,
            ty: wgpu::BindingType::Texture {
                sample_type,
                view_dimension: wgpu::TextureViewDimension::D2,
                multisampled: false,
            },
            count: None,
        };
        let uniform_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: None,
            entries: &[uniform_entry],
        });
        let post_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: None,
            entries: &[
                uniform_entry,
                texture_entry(1, wgpu::TextureSampleType::Float { filterable: false }),
                texture_entry(2, wgpu::TextureSampleType::Depth),
                texture_entry(3, wgpu::TextureSampleType::Float { filterable: false }),
            ],
        });
        let edge_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: None,
            entries: &[
                uniform_entry,
                texture_entry(2, wgpu::TextureSampleType::Depth),
            ],
        });
        let scene_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[Some(&uniform_layout)],
            immediate_size: 0,
        });
        let post_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[Some(&post_layout)],
            immediate_size: 0,
        });
        let edge_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[Some(&edge_layout)],
            immediate_size: 0,
        });

        let impostor_attributes = wgpu::vertex_attr_array![
            0 => Float32x4, 1 => Float32x4, 2 => Float32x4, 3 => Float32x4
        ];
        let impostor_buffers = [wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Impostor>() as u64,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &impostor_attributes,
        }];
        let vertex_attributes = wgpu::vertex_attr_array![
            0 => Float32x3, 1 => Float32x3, 2 => Float32x3
        ];
        let vertex_buffers = [wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &vertex_attributes,
        }];

        let depth = |write: bool| wgpu::DepthStencilState {
            format: DEPTH,
            depth_write_enabled: Some(write),
            depth_compare: Some(wgpu::CompareFunction::Less),
            stencil: Default::default(),
            bias: Default::default(),
        };
        let opaque_targets = [Some(wgpu::ColorTargetState::from(COLOR))];
        let pipeline = |vs: &str,
                        fs: &str,
                        buffers: &[wgpu::VertexBufferLayout],
                        layout: &wgpu::PipelineLayout,
                        targets: &[Option<wgpu::ColorTargetState>],
                        depth_stencil: Option<wgpu::DepthStencilState>| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some(fs),
                layout: Some(layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some(vs),
                    compilation_options: Default::default(),
                    buffers,
                },
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some(fs),
                    compilation_options: Default::default(),
                    targets,
                }),
                primitive: wgpu::PrimitiveState::default(),
                depth_stencil,
                multisample: Default::default(),
                multiview_mask: None,
                cache: None,
            })
        };

        let impostors = pipeline(
            "imp_vs",
            "imp_fs",
            &impostor_buffers,
            &scene_layout,
            &opaque_targets,
            Some(depth(true)),
        );
        let ribbons = pipeline(
            "tri_vs",
            "tri_fs",
            &vertex_buffers,
            &scene_layout,
            &opaque_targets,
            Some(depth(true)),
        );
        let glass = pipeline(
            "imp_vs",
            "glass_fs",
            &impostor_buffers,
            &scene_layout,
            &[Some(wgpu::ColorTargetState {
                format: OUT,
                blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                write_mask: wgpu::ColorWrites::ALL,
            })],
            Some(depth(false)),
        );
        let post = pipeline(
            "post_vs",
            "post_fs",
            &[],
            &post_pipeline_layout,
            &[Some(OUT.into())],
            None,
        );
        let edges = pipeline(
            "post_vs",
            "edge_fs",
            &[],
            &edge_pipeline_layout,
            &[Some(EDGE.into())],
            None,
        );

        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("uniforms"),
            size: std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Gpu {
            device,
            queue,
            max_size,
            uniform_layout,
            post_layout,
            edge_layout,
            impostors,
            glass,
            ribbons,
            edges,
            post,
            uniforms,
            state: Mutex::new(State::default()),
        }
    }

    /// Draw a frame, or `None` if it is larger than the device can draw.
    pub fn render(
        &self,
        scene: &Scene,
        cam: &Camera,
        width: usize,
        height: usize,
        page: [f32; 3],
        style: Style,
    ) -> Option<Frame> {
        let (w, h) = (width as u32, height as u32);
        if w == 0 || h == 0 || w > self.max_size || h > self.max_size {
            return None;
        }
        let mut state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        if state.scene != scene.version || scene.version == 0 {
            self.upload(&mut state, scene);
        }
        if state
            .targets
            .as_ref()
            .is_none_or(|t| t.width != w || t.height != h)
        {
            state.targets = Some(self.targets(w, h));
        }

        let m = cam.basis();
        let scale = (width.min(height) as f32) / (2.2 * cam.radius) * cam.zoom;
        let cx = width as f32 * 0.5 + cam.pan[0];
        let cy = height as f32 * 0.5 + cam.pan[1];

        // The depth the scene spans in this view, for the fog and for the
        // depth buffer's range; and the translucent layers in painter's order.
        let depth_of = |p: [f32; 3]| {
            m[2][0] * (p[0] - cam.center[0])
                + m[2][1] * (p[1] - cam.center[1])
                + m[2][2] * (p[2] - cam.center[2])
        };
        let (mut zmin, mut zmax) = (f32::MAX, f32::MIN);
        let mut track = |lo: f32, hi: f32| {
            zmin = zmin.min(lo);
            zmax = zmax.max(hi);
        };
        for s in &scene.spheres {
            if s.alpha > 0.0 && s.r > 0.0 {
                let z = depth_of(s.c);
                track(z - s.r, z + s.r);
            }
        }
        for c in &scene.capsules {
            if c.alpha > 0.0 {
                let (za, zb) = (depth_of(c.a), depth_of(c.b));
                track(za.min(zb) - c.r, za.max(zb) + c.r);
            }
        }
        for t in &scene.tris {
            for p in t.p {
                let z = depth_of(p);
                track(z, z);
            }
        }
        if zmin > zmax {
            (zmin, zmax) = (0.0, 0.0);
        }
        let span = (zmax - zmin).max(1e-3);

        let weight = (width.min(height) / 700).max(1);
        // Lines of one or two pixels are found by the ink pass itself.
        let fused = weight <= 2;
        let uniforms = Uniforms {
            rot: [
                [m[0][0], m[0][1], m[0][2], 0.0],
                [m[1][0], m[1][1], m[1][2], 0.0],
                [m[2][0], m[2][1], m[2][2], 0.0],
            ],
            center: [cam.center[0], cam.center[1], cam.center[2], 0.0],
            view: [cx, cy, scale, 0.0],
            size: [width as f32, height as f32, fused as u32 as f32, 0.0],
            depth: [zmax + 1.0, span + 2.0, zmin, span],
            page: [
                page[0],
                page[1],
                page[2],
                match style {
                    Style::Soft => 0.42,
                    Style::Illustrative => 0.0,
                },
            ],
            params: [
                (span * 0.012).max(0.45),
                match style {
                    Style::Soft => 0.0,
                    Style::Illustrative => 1.0,
                },
                weight as f32 + 0.5,
                weight as f32,
            ],
        };
        self.queue
            .write_buffer(&self.uniforms, 0, bytemuck::bytes_of(&uniforms));

        if !state.glass.is_empty() {
            // Sorted as integers: each layer's depth, its bits turned so they
            // order as the floats do, above its index, so a tie still sorts
            // the same way every frame.
            let mut order: Vec<u64> = state
                .glass
                .iter()
                .enumerate()
                .map(|(k, i)| {
                    let bits = ((depth_of(i.a) + depth_of(i.b)) * 0.5).to_bits();
                    let key = if bits >> 31 == 1 {
                        !bits
                    } else {
                        bits | 1 << 31
                    };
                    (key as u64) << 32 | k as u64
                })
                .collect();
            order.sort_unstable();
            let glass: Vec<Impostor> = order
                .iter()
                .map(|o| state.glass[(*o & 0xffff_ffff) as usize])
                .collect();
            let bytes: &[u8] = bytemuck::cast_slice(&glass);
            if state
                .glass_buffer
                .as_ref()
                .is_none_or(|(_, n)| *n < glass.len())
            {
                state.glass_buffer = Some((
                    self.device.create_buffer(&wgpu::BufferDescriptor {
                        label: Some("glass"),
                        size: bytes.len() as u64,
                        usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                        mapped_at_creation: false,
                    }),
                    glass.len(),
                ));
            }
            let (buffer, _) = state.glass_buffer.as_ref().unwrap();
            self.queue.write_buffer(buffer, 0, bytes);
        }

        let t = state.targets.as_ref().unwrap();
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let clear = |view, value| {
                Some(wgpu::RenderPassColorAttachment {
                    view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(value),
                        store: wgpu::StoreOp::Store,
                    },
                })
            };
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("scene"),
                color_attachments: &[clear(&t.color, wgpu::Color::TRANSPARENT)],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &t.depth,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_bind_group(0, &t.uniform, &[]);
            // Ribbons first: their depth is known before shading, so they
            // fill the depth buffer cheaply for the impostors to test against.
            if let Some((buffer, n)) = &state.ribbons {
                pass.set_pipeline(&self.ribbons);
                pass.set_vertex_buffer(0, buffer.slice(..));
                pass.draw(0..*n, 0..1);
            }
            if let Some((buffer, n)) = &state.opaque {
                pass.set_pipeline(&self.impostors);
                pass.set_vertex_buffer(0, buffer.slice(..));
                pass.draw(0..6, 0..*n);
            }
        }
        // Only heavy lines need their edges found in a pass of their own.
        if style == Style::Illustrative && !fused {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("edges"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &t.edge,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                ..Default::default()
            });
            pass.set_pipeline(&self.edges);
            pass.set_bind_group(0, &t.edges, &[]);
            pass.draw(0..3, 0..1);
        }
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("ink"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &t.out_view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                ..Default::default()
            });
            pass.set_pipeline(&self.post);
            pass.set_bind_group(0, &t.post, &[]);
            pass.draw(0..3, 0..1);
        }
        if let (false, Some((buffer, _))) = (state.glass.is_empty(), &state.glass_buffer) {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("glass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &t.out_view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &t.depth,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_pipeline(&self.glass);
            pass.set_bind_group(0, &t.uniform, &[]);
            pass.set_vertex_buffer(0, buffer.slice(..));
            pass.draw(0..6, 0..state.glass.len() as u32);
        }
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &t.out,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &t.readback,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(t.row),
                    rows_per_image: Some(h),
                },
            },
            wgpu::Extent3d {
                width: w,
                height: h,
                depth_or_array_layers: 1,
            },
        );
        let submission = self.queue.submit([encoder.finish()]);

        let slice = t.readback.slice(..);
        let mapped = std::sync::Arc::new(std::sync::OnceLock::new());
        let done = mapped.clone();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = done.set(r.is_ok());
        });
        // Asked to wait, Metal sleeps a millisecond at a time between checks,
        // which costs more than drawing the frame. Ask without waiting until
        // it is done.
        let start = std::time::Instant::now();
        while mapped.get().is_none() {
            self.device.poll(wgpu::PollType::Poll).ok()?;
            if start.elapsed() > std::time::Duration::from_secs(2) {
                return None;
            }
            std::thread::yield_now();
        }
        if mapped.get() != Some(&true) {
            return None;
        }
        let _ = submission;
        let mut bgra = vec![0u8; width * height * 4];
        {
            let mapped = slice.get_mapped_range();
            let line = width * 4;
            if t.row as usize == line {
                bgra.copy_from_slice(&mapped[..line * height]);
            } else {
                for (y, row) in bgra.chunks_exact_mut(line).enumerate() {
                    let o = y * t.row as usize;
                    row.copy_from_slice(&mapped[o..o + line]);
                }
            }
        }
        t.readback.unmap();
        Some(Frame {
            width,
            height,
            bgra,
        })
    }

    /// Put a new scene's primitives on the card.
    fn upload(&self, state: &mut State, scene: &Scene) {
        let mut opaque = Vec::new();
        let mut glass = Vec::new();
        for s in &scene.spheres {
            if s.alpha <= 0.0 || s.r <= 0.0 {
                continue;
            }
            let i = Impostor {
                a: s.c,
                r: s.r,
                b: s.c,
                kind: 0.0,
                ca: [s.color[0], s.color[1], s.color[2], s.alpha],
                cb: [s.color[0], s.color[1], s.color[2], s.alpha],
            };
            if s.alpha < OPAQUE {
                glass.push(i)
            } else {
                opaque.push(i)
            }
        }
        for c in &scene.capsules {
            if c.alpha <= 0.0 {
                continue;
            }
            let i = Impostor {
                a: c.a,
                r: c.r,
                b: c.b,
                kind: 1.0,
                ca: [c.ca[0], c.ca[1], c.ca[2], c.alpha],
                cb: [c.cb[0], c.cb[1], c.cb[2], c.alpha],
            };
            if c.alpha < OPAQUE {
                glass.push(i)
            } else {
                opaque.push(i)
            }
        }
        let vertices: Vec<Vertex> = scene
            .tris
            .iter()
            .flat_map(|t| {
                (0..3).map(move |k| Vertex {
                    p: t.p[k],
                    n: t.n[k],
                    c: t.c,
                })
            })
            .collect();

        let buffer = |label, bytes: &[u8]| {
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: bytes,
                    usage: wgpu::BufferUsages::VERTEX,
                })
        };
        state.opaque = (!opaque.is_empty()).then(|| {
            (
                buffer("impostors", bytemuck::cast_slice(&opaque)),
                opaque.len() as u32,
            )
        });
        state.ribbons = (!vertices.is_empty()).then(|| {
            (
                buffer("ribbons", bytemuck::cast_slice(&vertices)),
                vertices.len() as u32,
            )
        });
        state.glass = glass;
        state.scene = scene.version;
    }

    fn targets(&self, width: u32, height: u32) -> Targets {
        let texture = |label, format, usage| {
            self.device.create_texture(&wgpu::TextureDescriptor {
                label: Some(label),
                size: wgpu::Extent3d {
                    width,
                    height,
                    depth_or_array_layers: 1,
                },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage,
                view_formats: &[],
            })
        };
        let attachment = wgpu::TextureUsages::RENDER_ATTACHMENT;
        let read = attachment | wgpu::TextureUsages::TEXTURE_BINDING;
        let color = texture("color", COLOR, read).create_view(&Default::default());
        let depth = texture("depth", DEPTH, read).create_view(&Default::default());
        let edge = texture("edge", EDGE, read).create_view(&Default::default());
        let out = texture("out", OUT, attachment | wgpu::TextureUsages::COPY_SRC);
        let out_view = out.create_view(&Default::default());

        let row = (width * 4).div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT)
            * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: row as u64 * height as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let uniform = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.uniform_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: self.uniforms.as_entire_binding(),
            }],
        });
        let post = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.post_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.uniforms.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&color),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(&depth),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: wgpu::BindingResource::TextureView(&edge),
                },
            ],
        });
        let edges = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.edge_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.uniforms.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(&depth),
                },
            ],
        });
        Targets {
            width,
            height,
            color,
            depth,
            edge,
            out,
            out_view,
            readback,
            row,
            post,
            edges,
            uniform,
        }
    }
}
