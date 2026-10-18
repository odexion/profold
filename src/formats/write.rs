//! Writing a structure back out, as PDB or mmCIF: all of it, or the atoms
//! a filter keeps. Positions are those of the frame shown.

use std::fmt::Write;

use crate::protein::Structure;

/// The structure, or the atoms `keep` passes, as a PDB file.
pub fn pdb(s: &Structure, keep: impl Fn(usize) -> bool) -> String {
    let mut out = String::new();
    // Seventy characters a line, broken between words, as the reader joins
    // them with a space.
    let mut lines: Vec<String> = Vec::new();
    for word in s.title.split_whitespace() {
        match lines.last_mut() {
            Some(l) if l.len() + 1 + word.len() <= 70 => {
                l.push(' ');
                l.push_str(word);
            }
            _ => lines.push(word.to_string()),
        }
    }
    for (i, text) in lines.iter().enumerate() {
        if i == 0 {
            let _ = writeln!(out, "TITLE     {text}");
        } else {
            let _ = writeln!(out, "TITLE   {:>2}{text}", i + 1);
        }
    }
    let chains = pdb_chain_ids(s);
    let mut serial = 0usize;
    // In the order the atoms were read, which is the file's: a chain, its
    // TER, and the ligands after.
    let mut open: Option<(usize, String, i32)> = None;
    let ter = |out: &mut String, serial: &mut usize, (ci, name, seq): (usize, String, i32)| {
        *serial += 1;
        let _ = writeln!(out, "TER   {:>5}      {name:>3} {}{:>4}", *serial % 100_000, chains[ci], seq.rem_euclid(10_000));
    };
    for i in (0..s.atoms.len()).filter(|&i| keep(i)) {
        let a = &s.atoms[i];
        let res = &s.residues[a.residue as usize];
        let ci = res.chain as usize;
        let polymer = !res.outside_polymer();
        // A chain's polymer ends where another chain, or its ligands, begin.
        if let Some(prev) = open.take() {
            if prev.0 != ci || !polymer {
                ter(&mut out, &mut serial, prev);
            } else {
                open = Some(prev);
            }
        }
        serial += 1;
        let record = if polymer { "ATOM  " } else { "HETATM" };
        let symbol = if a.element.is_known() { a.element.symbol() } else { "" };
        let name = pdb_atom_name(a.name(), symbol);
        let res_name: String = res.name.chars().take(3).collect();
        let [x, y, z] = a.pos;
        let _ = writeln!(
            out,
            "{record}{:>5} {name} {res_name:>3} {}{:>4}    {x:>8.3}{y:>8.3}{z:>8.3}{:>6.2}{:>6.2}          {:>2}",
            serial % 100_000,
            chains[ci],
            res.seq.rem_euclid(10_000),
            1.0,
            a.b.clamp(-99.99, 999.99),
            symbol.to_uppercase(),
        );
        if polymer {
            open = Some((ci, res_name, res.seq));
        }
    }
    if let Some(prev) = open {
        ter(&mut out, &mut serial, prev);
    }
    out.push_str("END\n");
    out
}

/// PDB's atom name column: a one-letter element's names start a column in,
/// so the element lines up under a two-letter one's.
fn pdb_atom_name(name: &str, symbol: &str) -> String {
    let name: String = if name.is_empty() { symbol.to_uppercase() } else { name.chars().take(4).collect() };
    if name.len() < 4 && symbol.len() < 2 {
        format!(" {name:<3}")
    } else {
        format!("{name:<4}")
    }
}

/// One character per chain, as PDB has room for: a chain's own when it has
/// one to itself, else the next unused letter or digit.
fn pdb_chain_ids(s: &Structure) -> Vec<char> {
    const SPARE: &str = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let mut used: Vec<char> = Vec::new();
    let mut out = Vec::with_capacity(s.chains.len());
    for c in &s.chains {
        let own = c.id.chars().next().filter(|_| c.id.chars().count() == 1);
        let id = match own {
            Some(ch) if !used.contains(&ch) => ch,
            _ => SPARE.chars().find(|ch| !used.contains(ch) && !s.chains.iter().any(|o| o.id == ch.to_string())).unwrap_or(' '),
        };
        used.push(id);
        out.push(id);
    }
    out
}

/// The structure, or the atoms `keep` passes, as an mmCIF file.
pub fn mmcif(s: &Structure, keep: impl Fn(usize) -> bool) -> String {
    let mut out = String::new();
    let id: String = s.id.chars().filter(|c| c.is_ascii_alphanumeric() || *c == '_').collect();
    let _ = writeln!(out, "data_{}", if id.is_empty() { "profold" } else { &id });
    out.push_str("#\n");
    if !s.title.is_empty() {
        let _ = writeln!(out, "_struct.title {}", quote(&s.title));
        out.push_str("#\n");
    }
    write_assemblies(&mut out, s);
    out.push_str(
        "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
         _atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
         _atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n\
         _atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n\
         _atom_site.pdbx_PDB_model_num\n",
    );
    for (serial, i) in (0..s.atoms.len()).filter(|&i| keep(i)).enumerate() {
        let a = &s.atoms[i];
        let res = &s.residues[a.residue as usize];
        let chain = &s.chains[res.chain as usize];
        let asym = if chain.id.trim().is_empty() { "A" } else { chain.id.as_str() };
        let group = if res.outside_polymer() { "HETATM" } else { "ATOM" };
        let symbol = if a.element.is_known() { a.element.symbol() } else { "X" };
        let name = if a.name().is_empty() { symbol } else { a.name() };
        let [x, y, z] = a.pos;
        let _ = writeln!(
            out,
            "{group} {} {symbol} {} {} {} {} {x:.3} {y:.3} {z:.3} 1.00 {:.2} {} {} 1",
            serial + 1,
            quote(name),
            quote(&res.name),
            quote(asym),
            res.seq,
            a.b,
            res.seq,
            quote(asym),
        );
    }
    out.push_str("#\n");
    out
}

/// The assemblies, as `pdbx_struct_assembly_gen` and `pdbx_struct_oper_list`
/// loops, naming chains as the atom sites below do.
fn write_assemblies(out: &mut String, s: &Structure) {
    if s.assemblies.is_empty() {
        return;
    }
    let mut ops: Vec<crate::assembly::Op> = Vec::new();
    let mut rows: Vec<(String, String, String)> = Vec::new();
    for a in &s.assemblies {
        for part in &a.parts {
            let ids: Vec<String> = part
                .ops
                .iter()
                .map(|op| {
                    let i = ops.iter().position(|o| o == op).unwrap_or_else(|| {
                        ops.push(*op);
                        ops.len() - 1
                    });
                    (i + 1).to_string()
                })
                .collect();
            rows.push((a.id.clone(), ids.join(","), part.chains.join(",")));
        }
    }
    out.push_str("loop_\n_pdbx_struct_assembly.id\n_pdbx_struct_assembly.oligomeric_details\n");
    for a in &s.assemblies {
        let _ = writeln!(out, "{} {}", quote(&a.id), quote(&a.name));
    }
    out.push_str("#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n");
    for (id, expr, chains) in rows {
        let _ = writeln!(out, "{} {} {}", quote(&id), quote(&expr), quote(&chains));
    }
    out.push_str("#\nloop_\n_pdbx_struct_oper_list.id\n");
    for i in 1..=3 {
        for j in 1..=3 {
            let _ = writeln!(out, "_pdbx_struct_oper_list.matrix[{i}][{j}]");
        }
        let _ = writeln!(out, "_pdbx_struct_oper_list.vector[{i}]");
    }
    for (k, op) in ops.iter().enumerate() {
        let _ = write!(out, "{}", k + 1);
        for row in op {
            for v in row {
                let _ = write!(out, " {v:.6}");
            }
        }
        out.push('\n');
    }
    out.push_str("#\n");
}

/// A CIF value, quoted where it has to be.
fn quote(v: &str) -> String {
    if v.is_empty() {
        return ".".into();
    }
    let plain = !v.contains(char::is_whitespace)
        && !v.starts_with(['_', '#', '$', '\'', '"', ';', '[', ']'])
        && !matches!(v, "." | "?")
        && !v.to_ascii_lowercase().starts_with("data_")
        && !v.eq_ignore_ascii_case("loop_");
    if plain {
        v.to_string()
    } else if !v.contains('"') {
        format!("\"{v}\"")
    } else {
        format!("'{v}'")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::formats::{load_bytes, parse, Format};
    use crate::protein::Source;

    fn sample(name: &str) -> Structure {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples").join(name);
        load_bytes(name, std::fs::read(path).unwrap()).unwrap()
    }

    fn same(a: &Structure, b: &Structure) {
        assert_eq!(a.atoms.len(), b.atoms.len());
        assert_eq!(a.chains.len(), b.chains.len());
        assert_eq!(a.residues.len(), b.residues.len());
        for (x, y) in a.atoms.iter().zip(&b.atoms) {
            assert_eq!(x.name(), y.name());
            assert_eq!(x.element, y.element);
            assert!((0..3).all(|k| (x.pos[k] - y.pos[k]).abs() < 1e-3));
            assert!((x.b - y.b).abs() < 1e-2);
        }
        for (x, y) in a.residues.iter().zip(&b.residues) {
            assert_eq!((x.name.as_str(), x.seq, x.ca.is_some()), (y.name.as_str(), y.seq, y.ca.is_some()));
        }
    }

    #[test]
    fn pdb_and_mmcif_read_back_as_written() {
        let s = sample("4hhb.pdb.gz");
        let text = pdb(&s, |_| true);
        let back = parse(&text, Format::Pdb, "4hhb".into(), Source::File).unwrap();
        same(&s, &back);
        assert_eq!(back.title, s.title);
        let text = mmcif(&s, |_| true);
        let back = parse(&text, Format::Cif, "4hhb".into(), Source::File).unwrap();
        same(&s, &back);
        assert_eq!(back.title, s.title);
        assert_eq!(back.assemblies.len(), s.assemblies.len());
    }

    #[test]
    fn assemblies_are_written_with_the_atoms() {
        let mut s = sample("1crn.pdb.gz");
        let mut turn = crate::assembly::IDENTITY;
        turn[0][0] = -1.0;
        turn[1][1] = -1.0;
        turn[0][3] = 40.0;
        s.assemblies = vec![crate::assembly::Assembly {
            id: "1".into(),
            name: "dimeric".into(),
            parts: vec![crate::assembly::Part { chains: vec!["A".into()], ops: vec![crate::assembly::IDENTITY, turn] }],
        }];
        let back = parse(&mmcif(&s, |_| true), Format::Cif, "1crn".into(), Source::File).unwrap();
        assert_eq!(back.assemblies, s.assemblies);
    }

    #[test]
    fn a_filter_keeps_only_its_atoms() {
        let s = sample("1crn.pdb.gz");
        let text = pdb(&s, |i| s.atoms[i].residue < 5);
        let back = parse(&text, Format::Pdb, "part".into(), Source::File).unwrap();
        assert_eq!(back.residues.len(), 5);
        assert_eq!(back.atoms.len(), s.atoms.iter().filter(|a| a.residue < 5).count());
    }

    #[test]
    fn long_chain_names_get_letters_of_their_own() {
        let mut s = sample("1crn.pdb.gz");
        let mut copy = s.chains[0].clone();
        copy.id = "A2".into();
        s.chains.push(copy);
        assert_eq!(pdb_chain_ids(&s), vec!['A', 'B']);
        assert_eq!(pdb_atom_name("CA", "C"), " CA ");
        assert_eq!(pdb_atom_name("FE", "Fe"), "FE  ");
        assert_eq!(pdb_atom_name("HD21", "H"), "HD21");
        assert_eq!(quote("C1'"), "C1'");
        assert_eq!(quote("a title"), "\"a title\"");
    }
}
