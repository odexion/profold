//! Binary trajectories: CHARMM and NAMD's DCD, and GROMACS's compressed
//! XTC. Both hold positions alone, frame after frame, and are laid over a
//! structure that says what the atoms are.

use anyhow::{bail, Context, Result};

use super::cell_from_parameters;
use crate::protein::{Cell, Snapshot};

/// The atom count and every frame of a DCD file.
pub fn dcd(bytes: &[u8]) -> Result<(usize, Vec<Snapshot>)> {
    // Fortran records: a length, the data, the length again. The first
    // record is 84 bytes, which says which way round the numbers are.
    let little = match bytes.get(..4) {
        Some(b) if u32::from_le_bytes(b.try_into().unwrap()) == 84 => true,
        Some(b) if u32::from_be_bytes(b.try_into().unwrap()) == 84 => false,
        _ => bail!("not a DCD file"),
    };
    let int = |b: &[u8]| -> i32 {
        let a: [u8; 4] = b.try_into().unwrap();
        if little { i32::from_le_bytes(a) } else { i32::from_be_bytes(a) }
    };
    let mut at = 0usize;
    let mut record = || -> Option<&[u8]> {
        let len = int(bytes.get(at..at + 4)?) as usize;
        let body = bytes.get(at + 4..at + 4 + len)?;
        if bytes.get(at + 4 + len..at + 8 + len).map(int) != Some(len as i32) {
            return None;
        }
        at += len + 8;
        Some(body)
    };
    let head = record().context("the DCD header is truncated")?;
    if &head[..4] != b"CORD" {
        bail!("not a DCD file of coordinates");
    }
    let control: Vec<i32> = head[4..84].chunks_exact(4).map(int).collect();
    let fixed = control[8];
    let charmm = control[19] != 0;
    let has_cell = charmm && control[10] != 0;
    let four_d = charmm && control[11] != 0;
    if fixed > 0 {
        bail!("DCD files with fixed atoms are not supported");
    }
    record().context("the DCD title is missing")?;
    let natoms = int(record().context("the DCD atom count is missing")?.get(..4).context("no atom count")?) as usize;
    let floats = |b: &[u8]| -> Vec<f32> {
        b.chunks_exact(4)
            .map(|c| {
                let a: [u8; 4] = c.try_into().unwrap();
                if little { f32::from_le_bytes(a) } else { f32::from_be_bytes(a) }
            })
            .collect()
    };
    let mut frames = Vec::new();
    loop {
        let cell = if has_cell {
            let Some(b) = record() else { break };
            let v: Vec<f64> = b
                .chunks_exact(8)
                .map(|c| {
                    let a: [u8; 8] = c.try_into().unwrap();
                    if little { f64::from_le_bytes(a) } else { f64::from_be_bytes(a) }
                })
                .collect();
            dcd_cell(&v)
        } else {
            None
        };
        let (Some(x), Some(y), Some(z)) = (record(), record(), record()) else { break };
        let (x, y, z) = (floats(x), floats(y), floats(z));
        if x.len() != natoms || y.len() != natoms || z.len() != natoms {
            bail!("a DCD frame has the wrong number of atoms");
        }
        if four_d {
            record();
        }
        let pos = (0..natoms).map(|i| [x[i], y[i], z[i]]).collect();
        frames.push(Snapshot { step: None, pos, cell });
    }
    if frames.is_empty() {
        bail!("the DCD file has no frames");
    }
    Ok((natoms, frames))
}

/// CHARMM's cell record: A, γ, B, β, α, C, the angles as degrees or, from
/// newer writers, as their cosines.
fn dcd_cell(v: &[f64]) -> Option<Cell> {
    if v.len() < 6 || v[0] <= 0.0 || v[2] <= 0.0 || v[5] <= 0.0 {
        return None;
    }
    let angle = |x: f64| if x.abs() <= 1.0 { x.acos().to_degrees() } else { x } as f32;
    Some(cell_from_parameters(v[0] as f32, v[2] as f32, v[5] as f32, angle(v[4]), angle(v[3]), angle(v[1])))
}

/// The atom count and every frame of an XTC file, in ångström.
pub fn xtc(bytes: &[u8]) -> Result<(usize, Vec<Snapshot>)> {
    let mut r = Xdr { bytes, at: 0 };
    let mut frames = Vec::new();
    let mut natoms = 0;
    while r.at < bytes.len() {
        if r.int()? != 1995 {
            bail!("not an XTC file, or a frame in it is damaged");
        }
        let n = r.int()? as usize;
        if frames.is_empty() {
            natoms = n;
        } else if n != natoms {
            bail!("the atom count changes partway through the XTC file");
        }
        let _step = r.int()?;
        let _time = r.float()?;
        let b: Vec<f32> = (0..9).map(|_| r.float()).collect::<Result<_>>()?;
        let cell = (b[0] > 0.0).then(|| Cell {
            origin: [0.0; 3],
            vectors: [
                [b[0] * 10.0, b[1] * 10.0, b[2] * 10.0],
                [b[3] * 10.0, b[4] * 10.0, b[5] * 10.0],
                [b[6] * 10.0, b[7] * 10.0, b[8] * 10.0],
            ],
        });
        let coords = xtc_coords(&mut r, natoms)?;
        let pos = coords.chunks_exact(3).map(|c| [c[0] * 10.0, c[1] * 10.0, c[2] * 10.0]).collect();
        frames.push(Snapshot { step: None, pos, cell });
    }
    if frames.is_empty() {
        bail!("the XTC file has no frames");
    }
    Ok((natoms, frames))
}

/// XDR: big-endian, in four-byte words.
struct Xdr<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Xdr<'_> {
    fn word(&mut self) -> Result<[u8; 4]> {
        let w = self.bytes.get(self.at..self.at + 4).context("the XTC file is truncated")?;
        self.at += 4;
        Ok(w.try_into().unwrap())
    }

    fn int(&mut self) -> Result<i32> {
        Ok(i32::from_be_bytes(self.word()?))
    }

    fn float(&mut self) -> Result<f32> {
        Ok(f32::from_be_bytes(self.word()?))
    }

    fn opaque(&mut self, n: usize) -> Result<&[u8]> {
        let b = self.bytes.get(self.at..self.at + n).context("the XTC file is truncated")?;
        self.at += n.div_ceil(4) * 4;
        Ok(b)
    }
}

const MAGICINTS: [i32; 73] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64, 80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285,
    52015, 65536, 82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255, 1048576,
    1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216,
];
const FIRSTIDX: usize = 9;

/// The bits of the compressed coordinates, read most significant first, as
/// GROMACS's xdrfile writes them.
struct Bits<'a> {
    bytes: &'a [u8],
    at: usize,
    lastbits: u32,
    lastbyte: u32,
}

impl Bits<'_> {
    fn next_byte(&mut self) -> u32 {
        let b = self.bytes.get(self.at).copied().unwrap_or(0) as u32;
        self.at += 1;
        b
    }

    fn take(&mut self, mut n: u32) -> i32 {
        let mask = if n >= 32 { u32::MAX } else { (1u32 << n) - 1 };
        let mut num: u32 = 0;
        while n >= 8 {
            self.lastbyte = (self.lastbyte << 8) | self.next_byte();
            num |= (self.lastbyte >> self.lastbits) << (n - 8);
            n -= 8;
        }
        if n > 0 {
            if self.lastbits < n {
                self.lastbits += 8;
                self.lastbyte = (self.lastbyte << 8) | self.next_byte();
            }
            self.lastbits -= n;
            num |= (self.lastbyte >> self.lastbits) & ((1 << n) - 1);
        }
        (num & mask) as i32
    }

    /// Three numbers packed together into `bits` bits, each less than its size.
    fn ints(&mut self, bits: u32, sizes: [i32; 3]) -> [i32; 3] {
        let mut bytes = [0i32; 32];
        let mut count = 0;
        let mut bits = bits;
        while bits > 8 {
            bytes[count] = self.take(8);
            count += 1;
            bits -= 8;
        }
        if bits > 0 {
            bytes[count] = self.take(bits);
            count += 1;
        }
        let mut out = [0i32; 3];
        for i in (1..3).rev() {
            let mut num: i64 = 0;
            for j in (0..count).rev() {
                num = (num << 8) | bytes[j] as i64;
                let p = num / sizes[i] as i64;
                bytes[j] = p as i32;
                num -= p * sizes[i] as i64;
            }
            out[i] = num as i32;
        }
        out[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
        out
    }
}

fn sizeofint(size: i32) -> u32 {
    let (mut num, mut bits) = (1i64, 0u32);
    while size as i64 >= num && bits < 32 {
        bits += 1;
        num <<= 1;
    }
    bits
}

fn sizeofints(sizes: [i32; 3]) -> u32 {
    let mut bytes = [0u32; 32];
    bytes[0] = 1;
    let mut count = 1;
    for &size in &sizes {
        let mut tmp: u64 = 0;
        let mut k = 0;
        while k < count {
            tmp += bytes[k] as u64 * size as u64;
            bytes[k] = (tmp & 0xff) as u32;
            tmp >>= 8;
            k += 1;
        }
        while tmp != 0 {
            bytes[k] = (tmp & 0xff) as u32;
            k += 1;
            tmp >>= 8;
        }
        count = k;
    }
    let (mut num, mut bits) = (1u32, 0u32);
    count -= 1;
    while bytes[count] >= num {
        bits += 1;
        num *= 2;
    }
    bits + count as u32 * 8
}

/// One frame's coordinates, in nanometres, as xdrfile's
/// `xdrfile_decompress_coord_float` reads them.
fn xtc_coords(r: &mut Xdr, natoms: usize) -> Result<Vec<f32>> {
    let lsize = r.int()? as usize;
    if lsize != natoms {
        bail!("an XTC frame has the wrong number of atoms");
    }
    if lsize <= 9 {
        return (0..lsize * 3).map(|_| r.float()).collect();
    }
    let precision = r.float()?;
    let minint = [r.int()?, r.int()?, r.int()?];
    let maxint = [r.int()?, r.int()?, r.int()?];
    let sizeint: [i32; 3] = std::array::from_fn(|k| maxint[k] - minint[k] + 1);
    let (bitsize, bitsizeint) = if sizeint.iter().any(|&s| s > 0xffffff) {
        (0, sizeint.map(sizeofint))
    } else {
        (sizeofints(sizeint), [0; 3])
    };
    let mut smallidx = r.int()? as usize;
    if !(FIRSTIDX..MAGICINTS.len()).contains(&smallidx) {
        bail!("an XTC frame is damaged");
    }
    let mut smaller = MAGICINTS[FIRSTIDX.max(smallidx - 1)] / 2;
    let mut smallnum = MAGICINTS[smallidx] / 2;
    let mut sizesmall = [MAGICINTS[smallidx]; 3];
    let count = r.int()? as usize;
    let data = r.opaque(count)?;
    let mut bits = Bits { bytes: data, at: 0, lastbits: 0, lastbyte: 0 };
    let inv = 1.0 / precision;
    let mut out = Vec::with_capacity(lsize * 3);
    let mut run: i32 = 0;
    let mut i = 0;
    while i < lsize {
        let mut this = if bitsize == 0 {
            [bits.take(bitsizeint[0]), bits.take(bitsizeint[1]), bits.take(bitsizeint[2])]
        } else {
            bits.ints(bitsize, sizeint)
        };
        i += 1;
        for k in 0..3 {
            this[k] += minint[k];
        }
        let mut prev = this;
        let flag = bits.take(1);
        let mut is_smaller = 0;
        if flag == 1 {
            run = bits.take(5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller -= 1;
        }
        if run > 0 {
            for k in (0..run).step_by(3) {
                let mut small = bits.ints(smallidx as u32, sizesmall);
                i += 1;
                for j in 0..3 {
                    small[j] += prev[j] - smallnum;
                }
                if k == 0 {
                    // The first two of a run are written swapped, which packs
                    // a water's oxygen and hydrogens tighter.
                    std::mem::swap(&mut small, &mut prev);
                    out.extend(prev.iter().map(|&v| v as f32 * inv));
                } else {
                    prev = small;
                }
                out.extend(small.iter().map(|&v| v as f32 * inv));
            }
        } else {
            out.extend(this.iter().map(|&v| v as f32 * inv));
        }
        smallidx = (smallidx as i32 + is_smaller) as usize;
        if smallidx >= MAGICINTS.len() {
            bail!("an XTC frame is damaged");
        }
        if is_smaller < 0 {
            smallnum = smaller;
            smaller = if smallidx > FIRSTIDX { MAGICINTS[smallidx - 1] / 2 } else { 0 };
        } else if is_smaller > 0 {
            smaller = smallnum;
            smallnum = MAGICINTS[smallidx] / 2;
        }
        sizesmall = [MAGICINTS[smallidx]; 3];
    }
    out.truncate(lsize * 3);
    if out.len() != lsize * 3 {
        bail!("an XTC frame ends early");
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A DCD file as CHARMM writes one, little-endian, with a cell per frame.
    fn write_dcd(frames: &[Vec<[f32; 3]>]) -> Vec<u8> {
        let n = frames[0].len();
        let mut out = Vec::new();
        let mut record = |body: &[u8]| {
            out.extend((body.len() as i32).to_le_bytes());
            out.extend(body);
            out.extend((body.len() as i32).to_le_bytes());
        };
        let mut head = b"CORD".to_vec();
        let mut control = [0i32; 20];
        control[0] = frames.len() as i32;
        control[10] = 1;
        control[19] = 24;
        for c in control {
            head.extend(c.to_le_bytes());
        }
        record(&head);
        let mut title = 1i32.to_le_bytes().to_vec();
        title.extend([b' '; 80]);
        record(&title);
        record(&(n as i32).to_le_bytes());
        for f in frames {
            let cell: Vec<u8> = [20.0f64, 90.0, 30.0, 90.0, 90.0, 40.0].iter().flat_map(|v| v.to_le_bytes()).collect();
            record(&cell);
            for k in 0..3 {
                let axis: Vec<u8> = f.iter().flat_map(|p| p[k].to_le_bytes()).collect();
                record(&axis);
            }
        }
        out
    }

    #[test]
    fn a_dcd_reads_back_as_written() {
        let frames = vec![vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], vec![[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]]];
        let (n, read) = dcd(&write_dcd(&frames)).unwrap();
        assert_eq!(n, 2);
        assert_eq!(read.len(), 2);
        assert_eq!(read[1].pos, frames[1]);
        let cell = read[0].cell.unwrap();
        assert!((cell.vectors[0][0] - 20.0).abs() < 1e-4 && (cell.vectors[2][2] - 40.0).abs() < 1e-3);
        assert!(dcd(b"nothing of the sort").is_err());
    }

    #[test]
    fn bits_are_counted_as_xdrfile_counts_them() {
        assert_eq!(sizeofint(1), 1);
        assert_eq!(sizeofint(255), 8);
        assert_eq!(sizeofint(256), 9);
        // Bits enough to hold the product itself, one more than it needs:
        // what the files are written with, so what they are read with.
        assert_eq!(sizeofints([2, 2, 2]), 4);
        assert_eq!(sizeofints([256, 256, 256]), 25);
        let mut b = Bits { bytes: &[0b1011_0000, 0xff], at: 0, lastbits: 0, lastbyte: 0 };
        assert_eq!(b.take(1), 1);
        assert_eq!(b.take(3), 0b011);
        assert_eq!(b.take(8), 0b0000_1111);
    }

    /// MDAnalysis's adenylate kinase in water, downloaded into PROFOLD_TRAJ:
    /// the trajectory's first frame is the starting structure, to the XTC's
    /// precision.
    #[test]
    #[ignore]
    fn adk_xtc_starts_where_its_gro_does() {
        let Some(dir) = std::env::var_os("PROFOLD_TRAJ") else { return };
        let dir = std::path::Path::new(&dir);
        let (n, frames) = xtc(&std::fs::read(dir.join("adk_oplsaa.xtc")).unwrap()).unwrap();
        let gro = crate::formats::load_path(&dir.join("adk_oplsaa.gro")).unwrap();
        assert_eq!(n, gro.atoms.len());
        let worst = gro.atoms.iter().zip(&frames[0].pos).map(|(a, p)| (0..3).map(|k| (a.pos[k] - p[k]).abs()).fold(0.0, f32::max)).fold(0.0, f32::max);
        eprintln!("{} frames of {n} atoms; frame 0 within {worst} Å of the .gro", frames.len());
        assert!(worst < 0.011, "{worst}");
        let (dn, dframes) = dcd(&std::fs::read(dir.join("adk_dims.dcd")).unwrap()).unwrap();
        let pdb = crate::formats::load_path(&dir.join("adk_open.pdb")).unwrap();
        // Bonded neighbours stay bonded: the first atom's nearest is close.
        let p = &dframes[0].pos;
        let near = (1..dn).map(|i| ((p[0][0] - p[i][0]).powi(2) + (p[0][1] - p[i][1]).powi(2) + (p[0][2] - p[i][2]).powi(2)).sqrt()).fold(f32::MAX, f32::min);
        eprintln!("{} DCD frames of {dn} atoms (the PDB has {}); nearest neighbour {near} Å", dframes.len(), pdb.atoms.len());
        assert!(near > 0.8 && near < 1.6);
    }
}
