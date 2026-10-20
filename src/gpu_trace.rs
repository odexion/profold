//! The ray tracer's pixels on the graphics card, by a compute shader
//! (`trace.wgsl`) that follows the loop in `trace.rs`. The tree and the
//! primitives are built on the CPU and uploaded as they are; the picture is
//! traced in bands of rows, each read back before the next is sent, so no
//! one dispatch runs long enough to trip a driver's watchdog, the buffers
//! stay small, and the trace can be followed and given up between bands.

use std::sync::atomic::Ordering;

use wgpu::util::DeviceExt;

use super::Gpu;
use crate::trace::{Control, GpuNode, GpuParams, GpuPrim};

/// Samples traced in one dispatch.
const BAND: usize = 1 << 18;

pub(super) struct Tracer {
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
}

impl Tracer {
    pub(super) fn new(device: &wgpu::Device) -> Option<Tracer> {
        let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("profold trace"),
            source: wgpu::ShaderSource::Wgsl(include_str!("trace.wgsl").into()),
        });
        let storage = |binding, read_only| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("profold trace"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                storage(1, true),
                storage(2, true),
                storage(3, false),
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("profold trace"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("profold trace"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("trace"),
            compilation_options: Default::default(),
            cache: None,
        });
        if let Some(e) = pollster::block_on(scope.pop()) {
            eprintln!("profold: tracing on the CPU ({e})");
            return None;
        }
        Some(Tracer { pipeline, layout })
    }
}

impl Gpu {
    /// Trace `width` × `height` samples of the tree `nodes` over `prims`:
    /// each sample's shaded colour and depth, `f32::NEG_INFINITY` where
    /// nothing was met. `None` when the card cannot (the CPU traces
    /// instead), and `Some(None)` when `control` gives it up.
    pub fn trace(
        &self,
        nodes: &[GpuNode],
        prims: &[GpuPrim],
        params: GpuParams,
        width: usize,
        height: usize,
        control: &Control,
    ) -> Option<Option<(Vec<f32>, Vec<f32>)>> {
        let tracer = self
            .tracer
            .get_or_init(|| Tracer::new(&self.device))
            .as_ref()?;
        if nodes.is_empty() || prims.is_empty() || width == 0 || height == 0 {
            return None;
        }
        let rows = (BAND / width).clamp(1, height);
        let limit = self.device.limits().max_storage_buffer_binding_size as usize;
        let sizes = [
            std::mem::size_of_val(nodes),
            std::mem::size_of_val(prims),
            rows * width * 16,
        ];
        if sizes.iter().any(|&s| s > limit) {
            return None;
        }
        let device = &self.device;
        let usage = wgpu::BufferUsages::STORAGE;
        let nodes = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("trace nodes"),
            contents: bytemuck::cast_slice(nodes),
            usage,
        });
        let prims = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("trace prims"),
            contents: bytemuck::cast_slice(prims),
            usage,
        });
        let uniform = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trace params"),
            size: std::mem::size_of::<GpuParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let out_size = (rows * width * 16) as u64;
        let out = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trace out"),
            size: out_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trace readback"),
            size: out_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("trace"),
            layout: &tracer.layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: uniform.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: nodes.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: prims.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: out.as_entire_binding(),
                },
            ],
        });

        let mut color = vec![0.0f32; width * height * 3];
        let mut depth = vec![f32::NEG_INFINITY; width * height];
        let mut y0 = 0;
        while y0 < height {
            if control.cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
                return Some(None);
            }
            let n = rows.min(height - y0);
            let mut p = params;
            p.size = [width as u32, y0 as u32, n as u32, p.size[3]];
            self.queue.write_buffer(&uniform, 0, bytemuck::bytes_of(&p));
            let mut encoder = device.create_command_encoder(&Default::default());
            {
                let mut pass = encoder.begin_compute_pass(&Default::default());
                pass.set_pipeline(&tracer.pipeline);
                pass.set_bind_group(0, &bind, &[]);
                pass.dispatch_workgroups((width as u32).div_ceil(8), (n as u32).div_ceil(8), 1);
            }
            let bytes = (n * width * 16) as u64;
            encoder.copy_buffer_to_buffer(&out, 0, &readback, 0, bytes);
            self.queue.submit([encoder.finish()]);

            let slice = readback.slice(..bytes);
            let mapped = std::sync::Arc::new(std::sync::OnceLock::new());
            let done = mapped.clone();
            slice.map_async(wgpu::MapMode::Read, move |r| {
                let _ = done.set(r.is_ok());
            });
            let start = std::time::Instant::now();
            while mapped.get().is_none() {
                device.poll(wgpu::PollType::Poll).ok()?;
                if start.elapsed() > std::time::Duration::from_secs(30) {
                    return None;
                }
                std::thread::sleep(std::time::Duration::from_micros(200));
            }
            if mapped.get() != Some(&true) {
                return None;
            }
            {
                let data = slice.get_mapped_range();
                let samples: &[[f32; 4]] = bytemuck::cast_slice(&data);
                let at = y0 * width;
                for (k, s) in samples.iter().enumerate() {
                    if s[3] > -1.0e38 {
                        color[(at + k) * 3..(at + k) * 3 + 3].copy_from_slice(&s[..3]);
                        depth[at + k] = s[3];
                    }
                }
            }
            readback.unmap();
            if let Some(d) = control.done {
                d.fetch_add(n, Ordering::Relaxed);
            }
            y0 += n;
        }
        Some(Some((color, depth)))
    }
}
