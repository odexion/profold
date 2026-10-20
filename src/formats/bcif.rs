//! BinaryCIF, the PDB's and AlphaFold DB's compact mmCIF: categories of
//! columns, each column packed by a chain of encodings, all in MessagePack.
//!
//! The columns wanted are unpacked and written out as mmCIF text, which the
//! mmCIF reader then reads like any other.

use std::fmt::Write;

use anyhow::{anyhow, bail, Result};

use super::write::quote;

/// A MessagePack value, as much of it as BinaryCIF uses.
#[derive(Clone, Debug, PartialEq)]
enum Value {
    Nil,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    Bin(Vec<u8>),
    Array(Vec<Value>),
    Map(Vec<(Value, Value)>),
}

impl Value {
    fn get(&self, key: &str) -> Option<&Value> {
        match self {
            Value::Map(m) => m
                .iter()
                .find(|(k, _)| matches!(k, Value::Str(s) if s == key))
                .map(|(_, v)| v),
            _ => None,
        }
    }

    fn str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    fn int(&self) -> Option<i64> {
        match self {
            Value::Int(i) => Some(*i),
            Value::Float(f) => Some(*f as i64),
            _ => None,
        }
    }

    fn float(&self) -> Option<f64> {
        match self {
            Value::Int(i) => Some(*i as f64),
            Value::Float(f) => Some(*f),
            _ => None,
        }
    }

    fn array(&self) -> &[Value] {
        match self {
            Value::Array(a) => a,
            _ => &[],
        }
    }
}

struct Reader<'a> {
    b: &'a [u8],
    at: usize,
}

impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let s = self
            .b
            .get(self.at..self.at + n)
            .ok_or_else(|| anyhow!("the BinaryCIF file is truncated"))?;
        self.at += n;
        Ok(s)
    }

    fn uint(&mut self, n: usize) -> Result<u64> {
        Ok(self
            .take(n)?
            .iter()
            .fold(0u64, |acc, &b| (acc << 8) | b as u64))
    }

    fn value(&mut self) -> Result<Value> {
        let t = self.take(1)?[0];
        Ok(match t {
            0x00..=0x7f => Value::Int(t as i64),
            0x80..=0x8f => self.map((t & 0x0f) as usize)?,
            0x90..=0x9f => self.array((t & 0x0f) as usize)?,
            0xa0..=0xbf => self.string((t & 0x1f) as usize)?,
            0xc0 => Value::Nil,
            0xc2 => Value::Bool(false),
            0xc3 => Value::Bool(true),
            0xc4 => {
                let n = self.uint(1)? as usize;
                Value::Bin(self.take(n)?.to_vec())
            }
            0xc5 => {
                let n = self.uint(2)? as usize;
                Value::Bin(self.take(n)?.to_vec())
            }
            0xc6 => {
                let n = self.uint(4)? as usize;
                Value::Bin(self.take(n)?.to_vec())
            }
            0xca => Value::Float(f32::from_bits(self.uint(4)? as u32) as f64),
            0xcb => Value::Float(f64::from_bits(self.uint(8)?)),
            0xcc => Value::Int(self.uint(1)? as i64),
            0xcd => Value::Int(self.uint(2)? as i64),
            0xce => Value::Int(self.uint(4)? as i64),
            0xcf => Value::Int(self.uint(8)? as i64),
            0xd0 => Value::Int(self.uint(1)? as u8 as i8 as i64),
            0xd1 => Value::Int(self.uint(2)? as u16 as i16 as i64),
            0xd2 => Value::Int(self.uint(4)? as u32 as i32 as i64),
            0xd3 => Value::Int(self.uint(8)? as i64),
            0xd9 => {
                let n = self.uint(1)? as usize;
                self.string(n)?
            }
            0xda => {
                let n = self.uint(2)? as usize;
                self.string(n)?
            }
            0xdb => {
                let n = self.uint(4)? as usize;
                self.string(n)?
            }
            0xdc => {
                let n = self.uint(2)? as usize;
                self.array(n)?
            }
            0xdd => {
                let n = self.uint(4)? as usize;
                self.array(n)?
            }
            0xde => {
                let n = self.uint(2)? as usize;
                self.map(n)?
            }
            0xdf => {
                let n = self.uint(4)? as usize;
                self.map(n)?
            }
            0xe0..=0xff => Value::Int(t as i8 as i64),
            other => bail!("MessagePack type {other:#x} is not one BinaryCIF uses"),
        })
    }

    fn string(&mut self, n: usize) -> Result<Value> {
        Ok(Value::Str(
            String::from_utf8_lossy(self.take(n)?).into_owned(),
        ))
    }

    fn array(&mut self, n: usize) -> Result<Value> {
        (0..n)
            .map(|_| self.value())
            .collect::<Result<_>>()
            .map(Value::Array)
    }

    fn map(&mut self, n: usize) -> Result<Value> {
        (0..n)
            .map(|_| Ok((self.value()?, self.value()?)))
            .collect::<Result<_>>()
            .map(Value::Map)
    }
}

/// A column, unpacked: numbers, or text.
#[derive(Clone, Debug, PartialEq)]
enum Column {
    Ints(Vec<i64>),
    Floats(Vec<f64>),
    Strings(Vec<String>),
}

impl Column {
    fn ints(self) -> Result<Vec<i64>> {
        match self {
            Column::Ints(v) => Ok(v),
            Column::Floats(v) => Ok(v.into_iter().map(|x| x as i64).collect()),
            Column::Strings(_) => bail!("expected numeric data, found text"),
        }
    }

    fn len(&self) -> usize {
        match self {
            Column::Ints(v) => v.len(),
            Column::Floats(v) => v.len(),
            Column::Strings(v) => v.len(),
        }
    }

    fn text(&self, i: usize) -> String {
        match self {
            Column::Ints(v) => v[i].to_string(),
            Column::Floats(v) => {
                let x = v[i];
                if x.fract() == 0.0 && x.abs() < 1e9 {
                    format!("{x:.1}")
                } else {
                    format!("{x}")
                }
            }
            Column::Strings(v) => v[i].clone(),
        }
    }
}

/// Undo a column's encodings, the last applied first.
fn decode(data: &Value) -> Result<Column> {
    let mut current: Option<Column> = None;
    let mut raw: Option<&[u8]> = None;
    if let Some(Value::Bin(b)) = data.get("data") {
        raw = Some(b);
    }
    for e in data
        .get("encoding")
        .map(Value::array)
        .unwrap_or(&[])
        .iter()
        .rev()
    {
        let kind = e.get("kind").and_then(Value::str).unwrap_or("");
        current = Some(match kind {
            "ByteArray" => {
                let bytes = raw
                    .take()
                    .ok_or_else(|| anyhow!("a ByteArray with no bytes"))?;
                byte_array(bytes, e.get("type").and_then(Value::int).unwrap_or(0))?
            }
            "FixedPoint" => {
                let factor = e.get("factor").and_then(Value::float).unwrap_or(1.0);
                Column::Floats(
                    current
                        .take()
                        .ok_or_else(|| anyhow!("FixedPoint encoding has no input data"))?
                        .ints()?
                        .into_iter()
                        .map(|v| v as f64 / factor)
                        .collect(),
                )
            }
            "IntervalQuantization" => {
                let (min, max) = (
                    e.get("min").and_then(Value::float).unwrap_or(0.0),
                    e.get("max").and_then(Value::float).unwrap_or(1.0),
                );
                let steps = e.get("numSteps").and_then(Value::int).unwrap_or(2).max(2) as f64;
                let delta = (max - min) / (steps - 1.0);
                Column::Floats(
                    current
                        .take()
                        .ok_or_else(|| anyhow!("IntervalQuantization encoding has no input data"))?
                        .ints()?
                        .into_iter()
                        .map(|v| min + delta * v as f64)
                        .collect(),
                )
            }
            "RunLength" => {
                let pairs = current
                    .take()
                    .ok_or_else(|| anyhow!("RunLength encoding has no input data"))?
                    .ints()?;
                let mut out = Vec::new();
                for p in pairs.chunks_exact(2) {
                    out.extend(std::iter::repeat_n(p[0], p[1].max(0) as usize));
                }
                Column::Ints(out)
            }
            "Delta" => {
                let mut v = e.get("origin").and_then(Value::int).unwrap_or(0);
                let deltas = current
                    .take()
                    .ok_or_else(|| anyhow!("Delta encoding has no input data"))?
                    .ints()?;
                Column::Ints(
                    deltas
                        .into_iter()
                        .map(|d| {
                            v += d;
                            v
                        })
                        .collect(),
                )
            }
            "IntegerPacking" => {
                let packed = current
                    .take()
                    .ok_or_else(|| anyhow!("IntegerPacking encoding has no input data"))?
                    .ints()?;
                let unsigned = matches!(e.get("isUnsigned"), Some(Value::Bool(true)));
                let (upper, lower) = match (
                    e.get("byteCount").and_then(Value::int).unwrap_or(1),
                    unsigned,
                ) {
                    (1, true) => (0xff, i64::MIN),
                    (1, false) => (0x7f, -0x80),
                    (_, true) => (0xffff, i64::MIN),
                    (_, false) => (0x7fff, -0x8000),
                };
                let mut out = Vec::new();
                let mut j = 0;
                while j < packed.len() {
                    let mut value = 0;
                    let mut t = packed[j];
                    while t == upper || t == lower {
                        value += t;
                        j += 1;
                        t = *packed
                            .get(j)
                            .ok_or_else(|| anyhow!("IntegerPacking exceeds the data length"))?;
                    }
                    out.push(value + t);
                    j += 1;
                }
                Column::Ints(out)
            }
            "StringArray" => {
                let text = e.get("stringData").and_then(Value::str).unwrap_or("");
                let offsets = decode(&Value::Map(vec![
                    (
                        Value::Str("data".into()),
                        e.get("offsets").cloned().unwrap_or(Value::Nil),
                    ),
                    (
                        Value::Str("encoding".into()),
                        e.get("offsetEncoding")
                            .cloned()
                            .unwrap_or(Value::Array(vec![])),
                    ),
                ]))?
                .ints()?;
                let indices = decode(&Value::Map(vec![
                    (
                        Value::Str("data".into()),
                        raw.take()
                            .map(|b| Value::Bin(b.to_vec()))
                            .unwrap_or(Value::Nil),
                    ),
                    (
                        Value::Str("encoding".into()),
                        e.get("dataEncoding")
                            .cloned()
                            .unwrap_or(Value::Array(vec![])),
                    ),
                ]))?
                .ints()?;
                let words: Vec<&str> = offsets
                    .windows(2)
                    .map(|w| text.get(w[0] as usize..w[1] as usize).unwrap_or(""))
                    .collect();
                Column::Strings(
                    indices
                        .into_iter()
                        .map(|i| {
                            if i < 0 {
                                String::new()
                            } else {
                                words.get(i as usize).copied().unwrap_or("").to_string()
                            }
                        })
                        .collect(),
                )
            }
            other => bail!("the BinaryCIF encoding {other} is not known"),
        });
    }
    current.ok_or_else(|| anyhow!("a column with no data"))
}

fn byte_array(b: &[u8], kind: i64) -> Result<Column> {
    let c = |n: usize| b.chunks_exact(n);
    Ok(match kind {
        1 => Column::Ints(b.iter().map(|&x| x as i8 as i64).collect()),
        2 => Column::Ints(
            c(2).map(|x| i16::from_le_bytes([x[0], x[1]]) as i64)
                .collect(),
        ),
        3 => Column::Ints(
            c(4).map(|x| i32::from_le_bytes(x.try_into().unwrap()) as i64)
                .collect(),
        ),
        4 => Column::Ints(b.iter().map(|&x| x as i64).collect()),
        5 => Column::Ints(
            c(2).map(|x| u16::from_le_bytes([x[0], x[1]]) as i64)
                .collect(),
        ),
        6 => Column::Ints(
            c(4).map(|x| u32::from_le_bytes(x.try_into().unwrap()) as i64)
                .collect(),
        ),
        32 => Column::Floats(
            c(4).map(|x| f32::from_le_bytes(x.try_into().unwrap()) as f64)
                .collect(),
        ),
        33 => Column::Floats(
            c(8).map(|x| f64::from_le_bytes(x.try_into().unwrap()))
                .collect(),
        ),
        other => bail!("the BinaryCIF byte type {other} is not known"),
    })
}

/// The categories the mmCIF reader takes, written out as mmCIF text.
pub fn to_mmcif(bytes: &[u8]) -> Result<String> {
    const WANTED: [&str; 5] = [
        "_atom_site",
        "_struct",
        "_pdbx_struct_assembly",
        "_pdbx_struct_assembly_gen",
        "_pdbx_struct_oper_list",
    ];
    let file = Reader { b: bytes, at: 0 }.value()?;
    let block = file
        .get("dataBlocks")
        .map(Value::array)
        .and_then(|b| b.first())
        .ok_or_else(|| anyhow!("the BinaryCIF file has no data"))?;
    let header = block
        .get("header")
        .and_then(Value::str)
        .unwrap_or("profold");
    let mut out = format!("data_{header}\n#\n");
    for cat in block.get("categories").map(Value::array).unwrap_or(&[]) {
        let Some(name) = cat.get("name").and_then(Value::str) else {
            continue;
        };
        if !WANTED.contains(&name) {
            continue;
        }
        let rows = cat.get("rowCount").and_then(Value::int).unwrap_or(0) as usize;
        let mut columns: Vec<(String, Column, Option<Vec<i64>>)> = Vec::new();
        for col in cat.get("columns").map(Value::array).unwrap_or(&[]) {
            let Some(cname) = col.get("name").and_then(Value::str) else {
                continue;
            };
            let Some(data) = col.get("data") else {
                continue;
            };
            let column = decode(data)?;
            let mask = match col.get("mask") {
                Some(m @ Value::Map(_)) => Some(decode(m)?.ints()?),
                _ => None,
            };
            if column.len() == rows {
                columns.push((cname.to_string(), column, mask));
            }
        }
        if columns.is_empty() {
            continue;
        }
        let cell = |(_, c, mask): &(String, Column, Option<Vec<i64>>), i: usize| match mask
            .as_ref()
            .map(|m| m[i])
        {
            Some(1) => ".".to_string(),
            Some(2) => "?".to_string(),
            _ => quote(&c.text(i)),
        };
        if rows == 1 {
            for c in &columns {
                let _ = writeln!(out, "{name}.{} {}", c.0, cell(c, 0));
            }
        } else {
            out.push_str("loop_\n");
            for c in &columns {
                let _ = writeln!(out, "{name}.{}", c.0);
            }
            for i in 0..rows {
                let line: Vec<String> = columns.iter().map(|c| cell(c, i)).collect();
                out.push_str(&line.join(" "));
                out.push('\n');
            }
        }
        out.push_str("#\n");
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodings_undo_as_they_were_done() {
        let enc = |kinds: Vec<Value>| Value::Array(kinds);
        let map = |pairs: Vec<(&str, Value)>| {
            Value::Map(
                pairs
                    .into_iter()
                    .map(|(k, v)| (Value::Str(k.into()), v))
                    .collect(),
            )
        };
        // 1, 2, 3, 100: as deltas, packed into bytes, as a byte array.
        let bytes: Vec<u8> = [1i8, 1, 1, 97].iter().map(|&b| b as u8).collect();
        let column = map(vec![
            ("data", Value::Bin(bytes)),
            (
                "encoding",
                enc(vec![
                    map(vec![
                        ("kind", Value::Str("Delta".into())),
                        ("origin", Value::Int(0)),
                    ]),
                    map(vec![
                        ("kind", Value::Str("IntegerPacking".into())),
                        ("byteCount", Value::Int(1)),
                        ("isUnsigned", Value::Bool(false)),
                    ]),
                    map(vec![
                        ("kind", Value::Str("ByteArray".into())),
                        ("type", Value::Int(1)),
                    ]),
                ]),
            ),
        ]);
        assert_eq!(decode(&column).unwrap(), Column::Ints(vec![1, 2, 3, 100]));
        // A run of five sevens, then 12.5 by fixed point.
        let runs: Vec<u8> = [7i32, 5].iter().flat_map(|v| v.to_le_bytes()).collect();
        let column = map(vec![
            ("data", Value::Bin(runs)),
            (
                "encoding",
                enc(vec![
                    map(vec![("kind", Value::Str("RunLength".into()))]),
                    map(vec![
                        ("kind", Value::Str("ByteArray".into())),
                        ("type", Value::Int(3)),
                    ]),
                ]),
            ),
        ]);
        assert_eq!(decode(&column).unwrap(), Column::Ints(vec![7; 5]));
        // Packing past the byte's range carries on into the next.
        let packed: Vec<u8> = [127i8, 127, 46].iter().map(|&b| b as u8).collect();
        let column = map(vec![
            ("data", Value::Bin(packed)),
            (
                "encoding",
                enc(vec![
                    map(vec![
                        ("kind", Value::Str("IntegerPacking".into())),
                        ("byteCount", Value::Int(1)),
                        ("isUnsigned", Value::Bool(false)),
                    ]),
                    map(vec![
                        ("kind", Value::Str("ByteArray".into())),
                        ("type", Value::Int(1)),
                    ]),
                ]),
            ),
        ]);
        assert_eq!(decode(&column).unwrap(), Column::Ints(vec![300]));
    }

    /// RCSB's BinaryCIF, downloaded into PROFOLD_TRAJ, reads as its mmCIF does.
    #[test]
    #[ignore]
    fn rcsb_bcif_matches_its_mmcif() {
        let Some(dir) = std::env::var_os("PROFOLD_TRAJ") else {
            return;
        };
        let dir = std::path::Path::new(&dir);
        for id in ["1crn", "4hhb"] {
            let b = crate::formats::load_path(&dir.join(format!("{id}.bcif"))).unwrap();
            let sample = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(format!("samples/{id}.pdb.gz"));
            let p = crate::formats::load_path(&sample).unwrap();
            eprintln!(
                "{id}: bcif {} atoms {} chains '{}'; pdb {} atoms",
                b.atoms.len(),
                b.chains.len(),
                b.title,
                p.atoms.len()
            );
            assert_eq!(b.atoms.len(), p.atoms.len());
            assert_eq!(b.title.to_lowercase(), p.title.to_lowercase());
            for (x, y) in b.atoms.iter().zip(&p.atoms) {
                assert!(
                    (0..3).all(|k| (x.pos[k] - y.pos[k]).abs() < 1e-3),
                    "{:?} {:?}",
                    x.pos,
                    y.pos
                );
                assert_eq!(x.name(), y.name());
            }
        }
    }
}
