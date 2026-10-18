//! CIF in both its dialects: mmCIF for macromolecules, with Cartesian atom
//! sites, and the crystallographers' small-molecule CIF, with fractional sites
//! in the asymmetric unit and the symmetry operations that fill the cell.

use anyhow::{bail, Result};

use super::{cell_from_parameters, number};
use crate::protein::{Builder, Element, Record, Snapshot, Source, Structure};

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    if text.contains("_atom_site.") {
        parse_mmcif(text, id, source)
    } else {
        parse_small(text, id, source)
    }
}

fn unquote(v: &str) -> String {
    let v = v.trim();
    if v.len() >= 2
        && ((v.starts_with('\'') && v.ends_with('\'')) || (v.starts_with('"') && v.ends_with('"')))
    {
        v[1..v.len() - 1].to_string()
    } else {
        v.to_string()
    }
}

/// Split a CIF data line into values, honouring single and double quotes.
fn cif_tokens(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        while i < bytes.len() && (bytes[i] == b' ' || bytes[i] == b'\t') {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        let quote = bytes[i];
        if quote == b'\'' || quote == b'"' {
            i += 1;
            let start = i;
            // A quote only closes the value when whitespace follows it, so
            // O5' stays one name.
            while i < bytes.len()
                && !(bytes[i] == quote
                    && bytes.get(i + 1).is_none_or(|c| *c == b' ' || *c == b'\t'))
            {
                i += 1;
            }
            out.push(&line[start..i.min(bytes.len())]);
            i += 1;
        } else {
            let start = i;
            while i < bytes.len() && bytes[i] != b' ' && bytes[i] != b'\t' {
                i += 1;
            }
            out.push(&line[start..i]);
        }
    }
    out
}

fn parse_mmcif(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut b = Builder::new("mmCIF");
    let mut title = String::new();
    // Assemblies name chains by label_asym_id; the structure, by the
    // author's chain.
    let mut label_to_auth: std::collections::HashMap<String, String> = Default::default();

    let mut lines = text.lines().peekable();
    while let Some(line) = lines.next() {
        let trimmed = line.trim();

        if let Some(rest) = trimmed.strip_prefix("_struct.title") {
            if !rest.is_empty() && !rest.starts_with(char::is_whitespace) {
                continue;
            }
            let rest = rest.trim();
            if !rest.is_empty() {
                title = unquote(rest);
            } else if let Some(next) = lines.next() {
                // The value sits on the following line, either quoted or as a
                // semicolon-delimited block.
                let next = next.trim_end();
                if let Some(first) = next.strip_prefix(';') {
                    title = first.trim().to_string();
                    for line in lines.by_ref() {
                        let l = line.trim_end();
                        if l.trim() == ";" {
                            break;
                        }
                        if !title.is_empty() {
                            title.push(' ');
                        }
                        title.push_str(l.trim());
                    }
                } else {
                    title = unquote(next.trim());
                }
            }
            continue;
        }

        if trimmed != "loop_" {
            continue;
        }

        // Collect this loop's column names.
        let mut cols: Vec<String> = Vec::new();
        while let Some(peek) = lines.peek() {
            let p = peek.trim();
            if p.starts_with('_') {
                cols.push(p.split_whitespace().next().unwrap_or("").to_string());
                lines.next();
            } else {
                break;
            }
        }
        if !cols.first().is_some_and(|c| c.starts_with("_atom_site.")) {
            continue;
        }

        let idx = |name: &str| cols.iter().position(|c| c == name);
        let i_symbol = idx("_atom_site.type_symbol");
        let i_atom = idx("_atom_site.label_atom_id");
        let i_alt = idx("_atom_site.label_alt_id");
        let i_comp = idx("_atom_site.auth_comp_id").or_else(|| idx("_atom_site.label_comp_id"));
        let i_asym = idx("_atom_site.auth_asym_id").or_else(|| idx("_atom_site.label_asym_id"));
        let i_label_asym = idx("_atom_site.label_asym_id");
        let i_seq = idx("_atom_site.auth_seq_id").or_else(|| idx("_atom_site.label_seq_id"));
        let i_ins = idx("_atom_site.pdbx_PDB_ins_code");
        let i_x = idx("_atom_site.Cartn_x");
        let i_y = idx("_atom_site.Cartn_y");
        let i_z = idx("_atom_site.Cartn_z");
        let i_b = idx("_atom_site.B_iso_or_equiv");
        let i_model = idx("_atom_site.pdbx_PDB_model_num");

        let (Some(i_x), Some(i_y), Some(i_z), Some(i_comp), Some(i_asym), Some(i_seq)) =
            (i_x, i_y, i_z, i_comp, i_asym, i_seq)
        else {
            continue;
        };

        // Models after the first become frames.
        let mut first_model: Option<String> = None;
        let mut current_model: Option<String> = None;
        let mut frame: Vec<[f32; 3]> = Vec::new();
        for row in lines.by_ref() {
            let t = row.trim();
            if t.is_empty() || t == "#" || t.starts_with('_') || t == "loop_" {
                break;
            }
            let f = cif_tokens(row);
            if f.len() < cols.len() {
                continue;
            }
            let get = |i: Option<usize>| i.and_then(|i| f.get(i).copied()).unwrap_or("");

            let res_name = f[i_comp];
            if res_name == "HOH" || res_name == "DOD" {
                continue;
            }
            let alt = get(i_alt);
            if !(alt.is_empty() || alt == "." || alt == "?" || alt == "A") {
                continue;
            }
            let (Ok(x), Ok(y), Ok(z)) = (
                f[i_x].parse::<f32>(),
                f[i_y].parse::<f32>(),
                f[i_z].parse::<f32>(),
            ) else {
                continue;
            };

            if let Some(model) = i_model.and_then(|i| f.get(i).copied()) {
                if first_model.is_none() {
                    first_model = Some(model.to_string());
                    current_model = Some(model.to_string());
                }
                if current_model.as_deref() != Some(model) {
                    if !frame.is_empty() {
                        b.frames.push(Snapshot {
                            pos: std::mem::take(&mut frame),
                            cell: None,
                        });
                    }
                    current_model = Some(model.to_string());
                }
                if first_model.as_deref() != Some(model) {
                    frame.push([x, y, z]);
                    continue;
                }
            }

            let atom_name = get(i_atom);
            let element = {
                let sym = get(i_symbol);
                if sym.is_empty() || sym == "." {
                    Element::from_label(atom_name)
                } else {
                    Element::from_symbol(sym)
                }
            };
            let ins = match get(i_ins) {
                "?" | "." => "",
                v => v,
            };
            if let Some(label) = i_label_asym.and_then(|i| f.get(i)) {
                if !label_to_auth.contains_key(*label) {
                    label_to_auth.insert(label.to_string(), f[i_asym].to_string());
                }
            }
            b.push(Record {
                chain: f[i_asym],
                res_name,
                res_seq: f[i_seq].parse::<i32>().unwrap_or(0),
                ins,
                name: atom_name,
                element,
                pos: [x, y, z],
                b: get(i_b).parse::<f32>().unwrap_or(0.0),
                kind: 0,
            });
        }
        if !frame.is_empty() {
            b.frames.push(Snapshot { pos: frame, cell: None });
        }
        break;
    }

    if b.len() == 0 {
        bail!("no atoms found in mmCIF file");
    }
    b.assemblies = mmcif_assemblies(text, &label_to_auth);
    Ok(b.finish(id, title, source))
}

/// The assemblies an mmCIF file describes, from its `pdbx_struct_assembly`,
/// `_gen` and `oper_list` categories, which may be loops or single rows.
/// Only those categories are read, so the atom sites are not tokenised again.
fn mmcif_assemblies(text: &str, label_to_auth: &std::collections::HashMap<String, String>) -> Vec<crate::assembly::Assembly> {
    use crate::assembly::{expand, Assembly, Part};
    const WANTED: [&str; 3] = ["_pdbx_struct_assembly.", "_pdbx_struct_assembly_gen.", "_pdbx_struct_oper_list."];
    // Categories are separated by lines of a lone #.
    let mut picked = String::new();
    let mut block = String::new();
    let mut keep = false;
    for line in text.lines() {
        if line.trim() == "#" {
            if keep {
                picked.push_str(&block);
            }
            block.clear();
            keep = false;
            continue;
        }
        let t = line.trim_start();
        if t.starts_with('_') && WANTED.iter().any(|w| t.starts_with(w)) {
            keep = true;
        }
        block.push_str(line);
        block.push('\n');
    }
    if keep {
        picked.push_str(&block);
    }
    if picked.is_empty() {
        return Vec::new();
    }
    let (items, tables) = read_blocks(&picked);
    // Each category as rows of (column, value), loop or not.
    let rows = |prefix: &str| -> Vec<std::collections::HashMap<String, String>> {
        if let Some(t) = tables.iter().find(|t| t.cols.first().is_some_and(|c| c.starts_with(prefix))) {
            return t
                .rows
                .iter()
                .map(|r| {
                    t.cols
                        .iter()
                        .zip(r)
                        .map(|(c, v)| (c[prefix.len()..].to_string(), unquote(v)))
                        .collect()
                })
                .collect();
        }
        let one: std::collections::HashMap<String, String> = items
            .iter()
            .filter(|(k, _)| k.starts_with(prefix))
            .map(|(k, v)| (k[prefix.len()..].to_string(), unquote(v)))
            .collect();
        if one.is_empty() { Vec::new() } else { vec![one] }
    };
    let mut ops = std::collections::HashMap::new();
    for r in rows("_pdbx_struct_oper_list.") {
        let m = |k: &str| r.get(k).and_then(|v| v.parse::<f32>().ok()).unwrap_or(0.0);
        let op: crate::assembly::Op = std::array::from_fn(|i| {
            [
                m(&format!("matrix[{}][1]", i + 1)),
                m(&format!("matrix[{}][2]", i + 1)),
                m(&format!("matrix[{}][3]", i + 1)),
                m(&format!("vector[{}]", i + 1)),
            ]
        });
        if let Some(id) = r.get("id") {
            ops.insert(id.clone(), op);
        }
    }
    let names: std::collections::HashMap<String, String> = rows("_pdbx_struct_assembly.")
        .into_iter()
        .filter_map(|r| Some((r.get("id")?.clone(), r.get("oligomeric_details").cloned().unwrap_or_default())))
        .collect();
    let mut out: Vec<Assembly> = Vec::new();
    for r in rows("_pdbx_struct_assembly_gen.") {
        let (Some(id), Some(expr), Some(asyms)) = (r.get("assembly_id"), r.get("oper_expression"), r.get("asym_id_list")) else {
            continue;
        };
        let mut chains: Vec<String> = Vec::new();
        for label in asyms.split(',').map(str::trim) {
            if let Some(auth) = label_to_auth.get(label) {
                if !chains.contains(auth) {
                    chains.push(auth.clone());
                }
            }
        }
        let part = Part { chains, ops: expand(expr, &ops) };
        if part.chains.is_empty() || part.ops.is_empty() {
            continue;
        }
        match out.iter_mut().find(|a| &a.id == id) {
            Some(a) => a.parts.push(part),
            None => out.push(Assembly {
                id: id.clone(),
                name: names.get(id).cloned().unwrap_or_default(),
                parts: vec![part],
            }),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Small-molecule CIF
// ---------------------------------------------------------------------------

/// Every value in the file, as a flat list of tokens, with `;` text fields
/// folded into one token each. Enough to read tags and loops.
fn tokens(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        if let Some(first) = line.strip_prefix(';') {
            let mut field = first.to_string();
            for l in lines.by_ref() {
                if l.starts_with(';') {
                    break;
                }
                field.push(' ');
                field.push_str(l.trim());
            }
            out.push(field.trim().to_string());
            continue;
        }
        let line = match line.find(" #") {
            Some(i) => &line[..i],
            None if line.trim_start().starts_with('#') => "",
            None => line,
        };
        out.extend(cif_tokens(line).into_iter().map(str::to_string));
    }
    out
}

struct Table {
    cols: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn col(&self, names: &[&str]) -> Option<usize> {
        names
            .iter()
            .find_map(|n| self.cols.iter().position(|c| c.eq_ignore_ascii_case(n)))
    }
}

/// Tags with single values, and the loops, of the first data block that has
/// atom sites.
fn read_blocks(text: &str) -> (Vec<(String, String)>, Vec<Table>) {
    let toks = tokens(text);
    let mut items = Vec::new();
    let mut tables: Vec<Table> = Vec::new();
    let mut i = 0;
    let mut started = false;
    while i < toks.len() {
        let t = &toks[i];
        if t.to_ascii_lowercase().starts_with("data_") {
            // Stop at the second block once the first has given us atoms.
            if started
                && tables
                    .iter()
                    .any(|t| t.cols.iter().any(|c| c.starts_with("_atom_site_fract")))
            {
                break;
            }
            started = true;
            i += 1;
        } else if t.eq_ignore_ascii_case("loop_") {
            i += 1;
            let mut cols = Vec::new();
            while i < toks.len() && toks[i].starts_with('_') {
                cols.push(toks[i].clone());
                i += 1;
            }
            let mut values = Vec::new();
            while i < toks.len()
                && !toks[i].starts_with('_')
                && !toks[i].eq_ignore_ascii_case("loop_")
                && !toks[i].to_ascii_lowercase().starts_with("data_")
            {
                values.push(toks[i].clone());
                i += 1;
            }
            if !cols.is_empty() {
                let rows = values.chunks_exact(cols.len()).map(|c| c.to_vec()).collect();
                tables.push(Table { cols, rows });
            }
        } else if t.starts_with('_') {
            if let Some(v) = toks.get(i + 1) {
                items.push((t.clone(), v.clone()));
            }
            i += 2;
        } else {
            i += 1;
        }
    }
    (items, tables)
}

/// One symmetry operation, as a 3x4 affine map on fractional coordinates.
type Op = [[f32; 4]; 3];

/// "-x+1/2, y, -z+1/2" and its variants.
fn parse_op(s: &str) -> Option<Op> {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 3 {
        return None;
    }
    let mut op = [[0.0f32; 4]; 3];
    for (row, part) in parts.iter().enumerate() {
        let expr: String = part.chars().filter(|c| !c.is_whitespace()).collect();
        let expr = expr.to_ascii_lowercase();
        let mut sign = 1.0f32;
        let mut num = String::new();
        let flush = |num: &mut String, sign: f32, op: &mut Op| -> Option<()> {
            if num.is_empty() {
                return Some(());
            }
            let v = match num.split_once('/') {
                Some((a, b)) => a.parse::<f32>().ok()? / b.parse::<f32>().ok()?,
                None => num.parse::<f32>().ok()?,
            };
            op[row][3] += sign * v;
            num.clear();
            Some(())
        };
        for c in expr.chars() {
            match c {
                '+' | '-' => {
                    flush(&mut num, sign, &mut op)?;
                    sign = if c == '-' { -1.0 } else { 1.0 };
                }
                'x' | 'y' | 'z' => {
                    let k = (c as u8 - b'x') as usize;
                    // "2x" would be odd, but costs nothing to allow.
                    let scale = if num.is_empty() { 1.0 } else { num.parse().ok()? };
                    num.clear();
                    op[row][k] += sign * scale;
                    sign = 1.0;
                }
                '0'..='9' | '.' | '/' => num.push(c),
                _ => return None,
            }
        }
        flush(&mut num, sign, &mut op)?;
    }
    Some(op)
}

fn parse_small(text: &str, id: String, source: Source) -> Result<Structure> {
    let (items, tables) = read_blocks(text);
    let item = |name: &str| {
        items
            .iter()
            .find(|(k, _)| k.eq_ignore_ascii_case(name))
            .map(|(_, v)| v.as_str())
    };
    let cell_param = |name: &str| item(name).and_then(number);
    let cell = match (
        cell_param("_cell_length_a"),
        cell_param("_cell_length_b"),
        cell_param("_cell_length_c"),
    ) {
        (Some(a), Some(b), Some(c)) => Some(cell_from_parameters(
            a,
            b,
            c,
            cell_param("_cell_angle_alpha").unwrap_or(90.0),
            cell_param("_cell_angle_beta").unwrap_or(90.0),
            cell_param("_cell_angle_gamma").unwrap_or(90.0),
        )),
        _ => None,
    };

    let Some(sites) = tables.iter().find(|t| {
        t.col(&["_atom_site_fract_x", "_atom_site_cartn_x"]).is_some()
    }) else {
        bail!("no atom sites in CIF file");
    };
    let fractional = sites.col(&["_atom_site_fract_x"]).is_some();
    let (Some(ix), Some(iy), Some(iz)) = (
        sites.col(&["_atom_site_fract_x", "_atom_site_cartn_x"]),
        sites.col(&["_atom_site_fract_y", "_atom_site_cartn_y"]),
        sites.col(&["_atom_site_fract_z", "_atom_site_cartn_z"]),
    ) else {
        bail!("atom sites without coordinates");
    };
    let i_symbol = sites.col(&["_atom_site_type_symbol"]);
    let i_label = sites.col(&["_atom_site_label"]);

    let ops: Vec<Op> = tables
        .iter()
        .find_map(|t| {
            let c = t.col(&[
                "_symmetry_equiv_pos_as_xyz",
                "_space_group_symop_operation_xyz",
                "_space_group_symop.operation_xyz",
            ])?;
            Some(t.rows.iter().filter_map(|r| parse_op(&r[c])).collect())
        })
        .filter(|ops: &Vec<Op>| !ops.is_empty())
        .unwrap_or_else(|| vec![[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]]);

    let mut b = Builder::new("CIF");
    b.cell = cell;
    // Fractional positions already placed, for weeding out the copies that
    // symmetry puts on special positions.
    let mut placed: Vec<[f32; 3]> = Vec::new();
    for row in &sites.rows {
        let (Some(x), Some(y), Some(z)) = (number(&row[ix]), number(&row[iy]), number(&row[iz]))
        else {
            continue;
        };
        let label = i_label.map(|i| row[i].as_str()).unwrap_or("");
        let element = i_symbol
            .map(|i| {
                // Type symbols carry charges: "O2-", "Fe3+".
                let sym: String = row[i].chars().take_while(|c| c.is_ascii_alphabetic()).collect();
                Element::from_symbol(&sym)
            })
            .filter(|e| e.is_known())
            .unwrap_or_else(|| Element::from_label(label));
        let kind = b.kind(element.symbol());

        let (Some(cell), true) = (cell, fractional) else {
            b.push(Record {
                name: label,
                element,
                pos: [x, y, z],
                kind,
                ..Record::default()
            });
            continue;
        };
        for op in &ops {
            let f = [x, y, z];
            let mut g: [f32; 3] = std::array::from_fn(|r| {
                op[r][0] * f[0] + op[r][1] * f[1] + op[r][2] * f[2] + op[r][3]
            });
            for v in &mut g {
                *v = v.rem_euclid(1.0);
                if *v > 0.9999 {
                    *v = 0.0;
                }
            }
            let duplicate = placed.iter().any(|p| {
                let d: [f32; 3] = std::array::from_fn(|k| {
                    let d = g[k] - p[k];
                    d - d.round()
                });
                let c = cell.to_cartesian(d);
                let o = cell.origin;
                let v = [c[0] - o[0], c[1] - o[1], c[2] - o[2]];
                v[0] * v[0] + v[1] * v[1] + v[2] * v[2] < 0.01
            });
            if duplicate {
                continue;
            }
            placed.push(g);
            b.push(Record {
                name: label,
                element,
                pos: cell.to_cartesian(g),
                kind,
                ..Record::default()
            });
        }
    }
    if b.len() == 0 {
        bail!("no atoms found in CIF file");
    }
    let title = item("_chemical_name_mineral")
        .or_else(|| item("_chemical_name_common"))
        .or_else(|| item("_chemical_name_systematic"))
        .or_else(|| item("_chemical_formula_sum"))
        .filter(|t| *t != "?" && *t != ".")
        .unwrap_or("")
        .to_string();
    Ok(b.finish(id, title, source))
}

#[cfg(test)]
mod tests {
    use super::parse_op;

    #[test]
    fn symmetry_operations() {
        let op = parse_op("-x+1/2, y, -z+1/2").unwrap();
        assert_eq!(op[0], [-1.0, 0.0, 0.0, 0.5]);
        assert_eq!(op[1], [0.0, 1.0, 0.0, 0.0]);
        assert_eq!(op[2], [0.0, 0.0, -1.0, 0.5]);
        let op = parse_op("1/2+X,x-y,Z+0.25").unwrap();
        assert_eq!(op[0], [1.0, 0.0, 0.0, 0.5]);
        assert_eq!(op[1], [1.0, -1.0, 0.0, 0.0]);
        assert_eq!(op[2], [0.0, 0.0, 1.0, 0.25]);
    }
}
