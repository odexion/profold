//! Selections: a small query language for picking atoms out by what they
//! are and where, after PyMOL's and VMD's.
//!
//! ```text
//! chain A+B            resi 10-20,45         resn HIS+TYR
//! name CA              elem Fe               index 1-100
//! protein  ligand  water  backbone  sidechain  hydrogen
//! helix  strand  coil  all  none
//! within 5 of ligand   byres (name CA)       not water and chain A
//! type 1+2            c_pe < -3            vx >= 0.5             z > 10 and b < 30
//! ```
//!
//! A comparison reads a per-atom property of the file by its name, or `x`,
//! `y`, `z` and `b` (the B-factor).
//!
//! `and` binds tighter than `or`; `not`, `within` and `byres` tighter
//! still; parentheses group. Words are matched without regard to case.

use anyhow::{bail, Result};

use crate::protein::{is_solvent, Element, Ss, Structure};

#[derive(Clone, Debug, PartialEq)]
enum Expr {
    All,
    None,
    Chain(Vec<String>),
    Resi(Vec<(i32, i32)>),
    Resn(Vec<String>),
    Name(Vec<String>),
    Elem(Vec<String>),
    /// A simulation's atom types, by the names `Structure::kinds` gives them.
    Type(Vec<String>),
    Index(Vec<(i32, i32)>),
    Protein,
    Ligand,
    Water,
    Backbone,
    Sidechain,
    Hydrogen,
    Ss(Ss),
    Not(Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Within(f32, Box<Expr>),
    Byres(Box<Expr>),
    /// A value per atom against a number.
    Compare(String, Op, f32),
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Op {
    Less,
    LessEq,
    More,
    MoreEq,
    Equal,
    NotEqual,
}

impl Op {
    fn parse(t: &str) -> Option<Op> {
        Some(match t {
            "<" => Op::Less,
            "<=" => Op::LessEq,
            ">" => Op::More,
            ">=" => Op::MoreEq,
            "=" | "==" => Op::Equal,
            "!=" => Op::NotEqual,
            _ => return None,
        })
    }

    fn holds(self, a: f32, b: f32) -> bool {
        match self {
            Op::Less => a < b,
            Op::LessEq => a <= b,
            Op::More => a > b,
            Op::MoreEq => a >= b,
            Op::Equal => a == b,
            Op::NotEqual => a != b,
        }
    }
}

/// A selection, parsed and ready to be matched against a structure.
#[derive(Clone, Debug, PartialEq)]
pub struct Selection {
    text: String,
    expr: Expr,
}

impl Selection {
    pub fn parse(text: &str) -> Result<Self> {
        let tokens = tokenize(text);
        if tokens.is_empty() {
            bail!("an empty selection");
        }
        let mut p = Parser {
            tokens: &tokens,
            at: 0,
        };
        let expr = p.or()?;
        if let Some(t) = p.peek() {
            bail!("unexpected \"{t}\"");
        }
        Ok(Selection {
            text: text.trim().to_string(),
            expr,
        })
    }

    pub fn text(&self) -> &str {
        &self.text
    }

    /// Which atoms of `s` it takes, one flag per atom.
    pub fn atoms(&self, s: &Structure) -> Vec<bool> {
        eval(&self.expr, s)
    }
}

fn tokenize(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut word = String::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if matches!(c, '<' | '>' | '=' | '!') {
            // An operator, one character or two: `<`, `<=`, `!=`, `!`.
            if !word.is_empty() {
                out.push(std::mem::take(&mut word));
            }
            let mut op = c.to_string();
            if chars.peek() == Some(&'=') {
                op.push(chars.next().unwrap());
            }
            out.push(op);
        } else if c == '(' || c == ')' {
            if !word.is_empty() {
                out.push(std::mem::take(&mut word));
            }
            out.push(c.to_string());
        } else if c.is_whitespace() {
            if !word.is_empty() {
                out.push(std::mem::take(&mut word));
            }
        } else {
            word.push(c);
        }
    }
    if !word.is_empty() {
        out.push(word);
    }
    out
}

struct Parser<'a> {
    tokens: &'a [String],
    at: usize,
}

impl Parser<'_> {
    fn peek(&self) -> Option<&str> {
        self.tokens.get(self.at).map(String::as_str)
    }

    fn next(&mut self) -> Option<&str> {
        let t = self.tokens.get(self.at).map(String::as_str);
        self.at += 1;
        t
    }

    fn word(&self) -> Option<String> {
        self.peek().map(str::to_ascii_lowercase)
    }

    fn or(&mut self) -> Result<Expr> {
        let mut left = self.and()?;
        while self.word().as_deref() == Some("or") {
            self.next();
            left = Expr::Or(Box::new(left), Box::new(self.and()?));
        }
        Ok(left)
    }

    fn and(&mut self) -> Result<Expr> {
        let mut left = self.unary()?;
        loop {
            match self.word().as_deref() {
                Some("and") => {
                    self.next();
                    left = Expr::And(Box::new(left), Box::new(self.unary()?));
                }
                // Two terms side by side are both wanted, as in `chain A name CA`.
                Some(w) if w != "or" && w != ")" => {
                    left = Expr::And(Box::new(left), Box::new(self.unary()?));
                }
                _ => return Ok(left),
            }
        }
    }

    fn unary(&mut self) -> Result<Expr> {
        let Some(w) = self.word() else {
            bail!("incomplete selection expression")
        };
        match w.as_str() {
            "not" | "!" => {
                self.next();
                Ok(Expr::Not(Box::new(self.unary()?)))
            }
            "byres" | "byresidue" => {
                self.next();
                Ok(Expr::Byres(Box::new(self.unary()?)))
            }
            "within" | "around" => {
                self.next();
                let d: f32 = match self.next().map(str::parse::<f32>) {
                    Some(Ok(d)) if d >= 0.0 => d,
                    _ => bail!("within requires a distance in Å, as in \"within 5 of ligand\""),
                };
                if self.word().as_deref() == Some("of") {
                    self.next();
                }
                Ok(Expr::Within(d, Box::new(self.unary()?)))
            }
            "(" => {
                self.next();
                let e = self.or()?;
                if self.next() != Some(")") {
                    bail!("unmatched parenthesis");
                }
                Ok(e)
            }
            _ if self
                .tokens
                .get(self.at + 1)
                .and_then(|t| Op::parse(t))
                .is_some() =>
            {
                let name = self.next().unwrap_or_default().to_string();
                let op = Op::parse(self.next().unwrap_or_default()).expect("checked above");
                match self.next().map(str::parse::<f32>) {
                    Some(Ok(v)) => Ok(Expr::Compare(name, op, v)),
                    _ => bail!("{name} is compared with a number, as in \"{name} < 3\""),
                }
            }
            _ => self.term(),
        }
    }

    fn term(&mut self) -> Result<Expr> {
        let w = self.word().unwrap_or_default();
        self.next();
        let mut list = |what: &str| -> Result<Vec<String>> {
            match self.next() {
                Some(v) if v != "(" && v != ")" => Ok(v
                    .split(['+', ','])
                    .filter(|x| !x.is_empty())
                    .map(str::to_string)
                    .collect()),
                _ => bail!("{what} expects a value"),
            }
        };
        Ok(match w.as_str() {
            "all" | "*" => Expr::All,
            "none" => Expr::None,
            "chain" | "c." | "segid" => Expr::Chain(list("chain")?),
            "resn" | "resname" | "r." => Expr::Resn(list("resn")?),
            "name" | "n." => Expr::Name(list("name")?),
            "elem" | "element" | "e." => Expr::Elem(list("elem")?),
            "type" | "types" => Expr::Type(list("type")?),
            "resi" | "resid" | "i." => Expr::Resi(ranges(&list("resi")?)?),
            "index" | "id" => Expr::Index(ranges(&list("index")?)?),
            "protein" | "polymer" => Expr::Protein,
            "ligand" | "ligands" | "hetatm" | "organic" => Expr::Ligand,
            "water" | "solvent" | "waters" => Expr::Water,
            "backbone" | "bb." => Expr::Backbone,
            "sidechain" | "sc." => Expr::Sidechain,
            "hydrogen" | "hydrogens" | "h." => Expr::Hydrogen,
            "helix" | "helices" => Expr::Ss(Ss::Helix),
            "strand" | "sheet" | "strands" => Expr::Ss(Ss::Sheet),
            "coil" | "loop" | "loops" => Expr::Ss(Ss::Coil),
            "" => bail!("incomplete selection expression"),
            other => bail!("\"{other}\" is not a known selection keyword"),
        })
    }
}

/// `10-20,45` as inclusive ranges; a negative number keeps its sign.
fn ranges(items: &[String]) -> Result<Vec<(i32, i32)>> {
    items
        .iter()
        .map(|item| {
            let split = item
                .char_indices()
                .skip(1)
                .find(|&(_, c)| c == '-' || c == ':')
                .map(|(i, _)| i);
            let (a, b) = match split {
                Some(i) => (&item[..i], &item[i + 1..]),
                None => (item.as_str(), item.as_str()),
            };
            match (a.parse::<i32>(), b.parse::<i32>()) {
                (Ok(a), Ok(b)) => Ok((a.min(b), a.max(b))),
                _ => bail!("\"{item}\" is not a number or a range"),
            }
        })
        .collect()
}

fn eval(e: &Expr, s: &Structure) -> Vec<bool> {
    let n = s.atoms.len();
    let each = |f: &dyn Fn(usize) -> bool| (0..n).map(f).collect::<Vec<bool>>();
    let res = |i: usize| &s.residues[s.atoms[i].residue as usize];
    let any = |list: &[String], v: &str| list.iter().any(|x| x.eq_ignore_ascii_case(v));
    match e {
        Expr::All => vec![true; n],
        Expr::None => vec![false; n],
        Expr::Chain(ids) => each(&|i| {
            let id = s
                .chains
                .get(res(i).chain as usize)
                .map(|c| c.id.as_str())
                .unwrap_or("");
            any(ids, id)
        }),
        Expr::Resn(names) => each(&|i| any(names, &res(i).name)),
        Expr::Name(names) => each(&|i| any(names, s.atoms[i].name())),
        Expr::Type(ts) => each(&|i| {
            s.kinds
                .get(s.atoms[i].kind as usize)
                .is_some_and(|k| any(ts, k))
        }),
        Expr::Elem(es) => each(&|i| {
            let e = s.atoms[i].element;
            e.is_known() && es.iter().any(|x| Element::from_symbol(x) == e)
        }),
        Expr::Resi(rs) => each(&|i| {
            let seq = res(i).seq;
            rs.iter().any(|&(a, b)| (a..=b).contains(&seq))
        }),
        Expr::Index(rs) => each(&|i| rs.iter().any(|&(a, b)| (a..=b).contains(&(i as i32 + 1)))),
        Expr::Protein => each(&|i| res(i).ca.is_some()),
        Expr::Ligand => each(&|i| res(i).outside_polymer() && !is_solvent(&res(i).name)),
        Expr::Water => each(&|i| is_solvent(&res(i).name)),
        Expr::Backbone => each(&|i| {
            let r = res(i);
            r.ca.is_some() && [r.n, r.ca, r.c, r.o].contains(&Some(i as u32))
        }),
        Expr::Sidechain => each(&|i| {
            let r = res(i);
            // As PyMOL's: everything past the backbone, alpha carbon not included.
            r.ca.is_some() && ![r.n, r.ca, r.c, r.o].contains(&Some(i as u32))
        }),
        Expr::Hydrogen => each(&|i| s.atoms[i].element == Element::H),
        Expr::Ss(ss) => each(&|i| res(i).ca.is_some() && res(i).ss == *ss),
        Expr::Not(a) => eval(a, s).into_iter().map(|x| !x).collect(),
        Expr::And(a, b) => eval(a, s)
            .into_iter()
            .zip(eval(b, s))
            .map(|(x, y)| x && y)
            .collect(),
        Expr::Or(a, b) => eval(a, s)
            .into_iter()
            .zip(eval(b, s))
            .map(|(x, y)| x || y)
            .collect(),
        Expr::Byres(a) => {
            let inner = eval(a, s);
            let mut taken = vec![false; s.residues.len()];
            for (i, &x) in inner.iter().enumerate() {
                if x {
                    taken[s.atoms[i].residue as usize] = true;
                }
            }
            each(&|i| taken[s.atoms[i].residue as usize])
        }
        Expr::Within(d, a) => within(s, &eval(a, s), *d),
        Expr::Compare(name, op, v) => {
            let value: Box<dyn Fn(usize) -> f32> = match name.to_ascii_lowercase().as_str() {
                "x" => Box::new(|i| s.atoms[i].pos[0]),
                "y" => Box::new(|i| s.atoms[i].pos[1]),
                "z" => Box::new(|i| s.atoms[i].pos[2]),
                "b" | "bfactor" | "beta" => Box::new(|i| s.atoms[i].b),
                _ => match s.prop(name) {
                    Some(k) => {
                        let values = &s.props[k].values;
                        Box::new(move |i| values.get(i).copied().unwrap_or(f32::NAN))
                    }
                    None => return vec![false; n],
                },
            };
            each(&|i| op.holds(value(i), *v))
        }
    }
}

/// Atoms within `d` ångström of any atom marked, those included, found
/// through a grid of cells `d` wide.
fn within(s: &Structure, marked: &[bool], d: f32) -> Vec<bool> {
    use std::collections::HashMap;
    let cell = d.max(1.0);
    let key = |p: [f32; 3]| [0, 1, 2].map(|k| (p[k] / cell).floor() as i32);
    let mut grid: HashMap<[i32; 3], Vec<usize>> = HashMap::new();
    for (i, a) in s.atoms.iter().enumerate().filter(|(i, _)| marked[*i]) {
        grid.entry(key(a.pos)).or_default().push(i);
    }
    let d2 = d * d;
    (0..s.atoms.len())
        .map(|i| {
            if marked[i] {
                return true;
            }
            let p = s.atoms[i].pos;
            let c = key(p);
            (-1..=1).any(|dx| {
                (-1..=1).any(|dy| {
                    (-1..=1).any(|dz| {
                        grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz])
                            .is_some_and(|list| {
                                list.iter().any(|&j| {
                                    let q = s.atoms[j].pos;
                                    (p[0] - q[0]).powi(2)
                                        + (p[1] - q[1]).powi(2)
                                        + (p[2] - q[2]).powi(2)
                                        <= d2
                                })
                            })
                    })
                })
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(name: &str) -> Structure {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("samples")
            .join(name);
        crate::formats::load_bytes(name, std::fs::read(path).unwrap()).unwrap()
    }

    fn count(s: &Structure, q: &str) -> usize {
        Selection::parse(q)
            .unwrap_or_else(|e| panic!("{q}: {e}"))
            .atoms(s)
            .into_iter()
            .filter(|x| *x)
            .count()
    }

    #[test]
    fn comparisons_read_properties_positions_and_b_factors() {
        let dump = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS pp pp pp\n0 4\n0 4\n0 4\n\
                    ITEM: ATOMS id type x y z c_pe\n1 1 0 0 0 -4.0\n2 1 1 1 1 -2.5\n3 2 2 2 3 -1.0\n";
        let s = crate::formats::load_bytes("t.lammpstrj", dump.as_bytes().to_vec()).unwrap();
        assert_eq!(count(&s, "c_pe < -3"), 1);
        assert_eq!(count(&s, "c_pe<=-2.5"), 2);
        assert_eq!(count(&s, "C_PE >= -2.5 and z > 2"), 1);
        assert_eq!(count(&s, "c_pe != -1"), 2);
        assert_eq!(count(&s, "! c_pe = -4"), 2);
        assert_eq!(count(&s, "nothing_like_it > 0"), 0);
        assert!(Selection::parse("c_pe < high").is_err());
        let crambin = sample("1crn.pdb.gz");
        let all = crambin.atoms.len();
        assert_eq!(count(&crambin, "b < 1000"), all);
        assert_eq!(count(&crambin, "x > 1e9"), 0);
    }

    #[test]
    fn terms_pick_what_they_name() {
        let s = sample("4hhb.pdb.gz");
        let n = s.atoms.len();
        assert_eq!(count(&s, "all"), n);
        assert_eq!(count(&s, "none"), 0);
        let chain_a = s
            .atoms
            .iter()
            .filter(|a| s.chains[s.residues[a.residue as usize].chain as usize].id == "A")
            .count();
        assert_eq!(count(&s, "chain A"), chain_a);
        assert_eq!(count(&s, "chain A+B"), count(&s, "chain A or chain B"));
        assert_eq!(
            count(&s, "name CA"),
            s.residues.iter().filter(|r| r.ca.is_some()).count()
        );
        assert_eq!(count(&s, "elem Fe"), 4);
        assert_eq!(count(&s, "resn HEM and elem Fe"), 4);
        assert_eq!(count(&s, "ligand and elem Fe"), 4);
        assert_eq!(count(&s, "not all"), 0);
        assert_eq!(
            count(&s, "protein") + count(&s, "ligand") + count(&s, "water"),
            n
        );
        assert_eq!(
            count(&s, "backbone") + count(&s, "sidechain"),
            count(&s, "protein")
        );
        assert_eq!(count(&s, "resi 1-10 and chain A and name CA"), 10);
        assert_eq!(
            count(&s, "chain A resi 1-10 name CA"),
            10,
            "side by side is and"
        );
        assert_eq!(count(&s, "index 1-5"), 5);
        assert!(count(&s, "helix") > count(&s, "strand"));
    }

    #[test]
    fn distance_and_whole_residues() {
        let s = sample("4hhb.pdb.gz");
        let near = count(&s, "within 5 of resn HEM");
        assert!(near > count(&s, "resn HEM"));
        // The proximal histidine is near its haem; byres takes all of it.
        assert_eq!(
            count(&s, "byres (resn HIS and chain A and resi 87)"),
            count(&s, "chain A and resi 87")
        );
        let his_near = count(&s, "byres (resn HIS and within 3 of elem Fe)");
        assert_eq!(
            his_near % 10,
            0,
            "whole histidines, ten atoms each: {his_near}"
        );
        assert!(his_near >= 40);
        assert_eq!(count(&s, "within 0 of elem Fe"), 4);
    }

    #[test]
    fn precedence_and_errors() {
        let s = sample("1crn.pdb.gz");
        assert_eq!(
            count(&s, "name CA or name N and resi 1"),
            count(&s, "name CA or (name N and resi 1)")
        );
        assert_eq!(
            count(&s, "not name CA and resi 1"),
            count(&s, "(not name CA) and resi 1")
        );
        for bad in [
            "",
            "chain",
            "resi x",
            "within of ligand",
            "(chain A",
            "frobnicate",
            "chain A )",
        ] {
            assert!(Selection::parse(bad).is_err(), "{bad:?} parsed");
        }
        assert_eq!(count(&s, "RESN thr"), count(&s, "resn THR"));
    }
}
