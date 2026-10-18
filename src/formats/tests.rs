//! One small file per format, read back and checked.

use super::*;

fn read(text: &str, format: Format) -> Structure {
    parse(text, format, "test".into(), Source::File).unwrap()
}

fn close(a: [f32; 3], b: [f32; 3]) -> bool {
    (0..3).all(|k| (a[k] - b[k]).abs() < 1e-3)
}

#[test]
fn names() {
    assert_eq!(Format::from_name("1crn.pdb.gz"), Some(Format::Pdb));
    assert_eq!(Format::from_name("data.melt"), Some(Format::LammpsData));
    assert_eq!(Format::from_name("dump.crack"), Some(Format::LammpsDump));
    assert_eq!(Format::from_name("traj.lammpstrj"), Some(Format::LammpsDump));
    assert_eq!(Format::from_name("CONTCAR"), Some(Format::Poscar));
    assert_eq!(Format::from_name("XDATCAR"), Some(Format::Poscar));
    assert_eq!(Format::from_name("conf.gro"), Some(Format::Gro));
    assert_eq!(Format::from_name("notes.txt"), None);
}

const WATER_XYZ: &str = "3
water
O 0.000 0.000 0.000
H 0.957 0.000 0.000
H -0.240 0.927 0.000
3
water, later
O 0.000 0.000 0.100
H 0.957 0.000 0.100
H -0.240 0.927 0.100
";

#[test]
fn xyz_frames_and_bonds() {
    let s = read(WATER_XYZ, Format::Xyz);
    assert_eq!(s.atoms.len(), 3);
    assert_eq!(s.title, "water");
    assert_eq!(s.frame_count(), 2);
    assert_eq!(s.bonds.len(), 2);
    assert_eq!(s.atoms[0].element, Element::O);
    assert_eq!(Format::sniff(WATER_XYZ), Some(Format::Xyz));
}

#[test]
fn extended_xyz() {
    let text = "2
Lattice=\"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"
Cu 0.0 0.0 0.0
Cu 2.0 2.0 0.0
";
    let s = read(text, Format::Xyz);
    assert_eq!(s.cell.unwrap().lengths(), [4.0, 4.0, 4.0]);
    assert_eq!(s.atoms[1].element, Element::from_symbol("Cu"));
    // Copper does not bond to copper.
    assert!(s.bonds.is_empty());
}

#[test]
fn gro() {
    let text = "Two waters t= 0.0
    6
    1SOL     OW    1   0.126   0.163   0.144
    1SOL    HW1    2   0.190   0.227   0.183
    1SOL    HW2    3   0.177   0.096   0.098
    2NA      NA    4   1.000   1.000   1.000
    3CL      CL    5   0.500   0.500   0.500
    4SOL     OW    6   0.300   0.300   0.300
   2.00000   2.00000   2.00000
";
    let s = read(text, Format::Gro);
    assert_eq!(s.title, "Two waters");
    assert_eq!(s.atoms.len(), 6);
    assert!(close(s.atoms[0].pos, [1.26, 1.63, 1.44]));
    assert_eq!(s.atoms[3].element, Element::from_symbol("Na"));
    assert_eq!(s.atoms[4].element, Element::from_symbol("Cl"));
    assert_eq!(s.cell.unwrap().lengths(), [20.0, 20.0, 20.0]);
    assert_eq!(s.bonds.len(), 2);
}

#[test]
fn pdb_models_become_frames() {
    let text = "\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1       0.100   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.558   0.000   0.000  1.00  0.00           C
ENDMDL
END
";
    let s = read(text, Format::Pdb);
    assert_eq!(s.atoms.len(), 2);
    assert_eq!(s.frame_count(), 2);
    assert!(s.cell.is_some());
    assert_eq!(s.bonds, vec![[0, 1]]);
}

#[test]
fn pqr() {
    let text = "\
ATOM      1  N   ALA     1      -0.677  -1.230  -0.491 -0.3000 1.8240
ATOM      2  CA  ALA A   1      -0.001   0.064  -0.491  0.0337 1.9080
";
    let s = read(text, Format::Pqr);
    assert_eq!(s.atoms.len(), 2);
    assert!(close(s.atoms[1].pos, [-0.001, 0.064, -0.491]));
}

#[test]
fn mol2() {
    let text = "\
@<TRIPOS>MOLECULE
methanol
 6 5 1 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3     1  MOH1        0.0000
      2 O1          1.4300    0.0000    0.0000 O.3     1  MOH1        0.0000
      3 H1         -0.3600    1.0300    0.0000 H       1  MOH1        0.0000
      4 H2         -0.3600   -0.5100    0.8900 H       1  MOH1        0.0000
      5 H3         -0.3600   -0.5100   -0.8900 H       1  MOH1        0.0000
      6 H4          1.7500    0.9100    0.0000 H       1  MOH1        0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
     3     1     4    1
     4     1     5    1
     5     2     6    1
";
    let s = read(text, Format::Mol2);
    assert_eq!(s.title, "methanol");
    assert_eq!(s.atoms.len(), 6);
    assert_eq!(s.bonds.len(), 5);
    assert_eq!(s.atoms[1].element, Element::O);
    assert_eq!(s.residues[0].name, "MOH");
}

#[test]
fn sdf_v2000() {
    let text = "\
water
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9570    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
$$$$
";
    let s = read(text, Format::Sdf);
    assert_eq!(s.title, "water");
    assert_eq!(s.atoms.len(), 3);
    assert_eq!(s.bonds, vec![[0, 1], [0, 2]]);
    assert_eq!(Format::sniff(text), Some(Format::Sdf));
}

#[test]
fn lammps_data() {
    let text = "LAMMPS data file

3 atoms
1 bonds
2 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 15.9994
2 1.008

Atoms # full

3 1 2 0.41 0.957 0.0 0.0
1 1 1 -0.82 0.0 0.0 0.0
2 1 2 0.41 -0.24 0.927 0.0

Bonds

1 1 1 3
";
    let s = read(text, Format::LammpsData);
    assert_eq!(s.atoms.len(), 3);
    assert_eq!(s.atoms[0].element, Element::O);
    assert_eq!(s.atoms[2].element, Element::H);
    assert!(close(s.atoms[2].pos, [0.957, 0.0, 0.0]));
    assert_eq!(s.bonds, vec![[0, 2]]);
    assert_eq!(s.cell.unwrap().lengths(), [10.0, 10.0, 10.0]);
    assert_eq!(Format::sniff(text), Some(Format::LammpsData));
}

#[test]
fn lammps_dump() {
    let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 2.0
0.0 2.0
0.0 2.0
ITEM: ATOMS id type xs ys zs
2 2 0.5 0.5 0.5
1 1 0.0 0.0 0.0
ITEM: TIMESTEP
100
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 2.0
0.0 2.0
0.0 2.0
ITEM: ATOMS id type xs ys zs
1 1 0.1 0.0 0.0
2 2 0.5 0.5 0.5
";
    let s = read(text, Format::LammpsDump);
    assert_eq!(s.atoms.len(), 2);
    assert_eq!(s.frame_count(), 2);
    assert!(close(s.atoms[1].pos, [1.0, 1.0, 1.0]));
    assert!(close(s.frames[1].pos[0], [0.2, 0.0, 0.0]));
    assert_eq!(s.kinds, vec!["1".to_string(), "2".to_string()]);
    assert!(!s.atoms[0].element.is_known());
    // Atoms of no element are drawn half their spacing across.
    assert!(s.unknown_radius < 1.6);
}

#[test]
fn poscar() {
    let text = "Si diamond
1.0
5.43 0.0 0.0
0.0 5.43 0.0
0.0 0.0 5.43
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
";
    let s = read(text, Format::Poscar);
    assert_eq!(s.atoms.len(), 2);
    assert!(close(s.atoms[1].pos, [1.3575, 1.3575, 1.3575]));
    assert_eq!(s.atoms[0].element, Element::from_symbol("Si"));
    assert_eq!(s.bonds.len(), 1);
}

#[test]
fn small_molecule_cif_fills_the_cell() {
    // Rock salt: two sites and the face-centring operations.
    let text = "data_NaCl
_cell_length_a 5.64
_cell_length_b 5.64
_cell_length_c 5.64
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_chemical_name_mineral 'Halite'
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
x,y+1/2,z+1/2
x+1/2,y,z+1/2
x+1/2,y+1/2,z
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1 Na+ 0.0 0.0 0.0
Cl1 Cl- 0.5 0.5 0.5
";
    let s = read(text, Format::Cif);
    assert_eq!(s.title, "Halite");
    assert_eq!(s.atoms.len(), 8);
    assert_eq!(
        s.atoms.iter().filter(|a| a.element == Element::from_symbol("Na")).count(),
        4
    );
}

#[test]
fn cube() {
    let text = "title
comment
    2    0.000000    0.000000    0.000000
    2    1.000000    0.000000    0.000000
    2    0.000000    1.000000    0.000000
    2    0.000000    0.000000    1.000000
    8    8.000000    0.000000    0.000000    0.000000
    1    1.000000    1.808846    0.000000    0.000000
";
    let s = read(text, Format::Cube);
    assert_eq!(s.atoms.len(), 2);
    assert!(close(s.atoms[1].pos, [0.9572, 0.0, 0.0]));
}

#[test]
fn xsf() {
    let text = "CRYSTAL
PRIMVEC
 2.0 0.0 0.0
 0.0 2.0 0.0
 0.0 0.0 2.0
PRIMCOORD
 2 1
 Na 0.0 0.0 0.0
 17 1.0 1.0 1.0
";
    let s = read(text, Format::Xsf);
    assert_eq!(s.atoms.len(), 2);
    assert_eq!(s.atoms[1].element, Element::from_symbol("Cl"));
    assert!(s.cell.is_some());
}

#[test]
fn gzip() {
    use std::io::Write;
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    enc.write_all(WATER_XYZ.as_bytes()).unwrap();
    let text = decode(enc.finish().unwrap()).unwrap();
    assert_eq!(text, WATER_XYZ);
}

const SER: &str = "\
ATOM      1  N   SER A  42      11.104   6.134  -6.504  1.00 12.50           N
ATOM      2  CA  SER A  42      11.639   6.071  -5.147  1.00 13.25           C
ATOM      3  C   SER A  42       9.766   6.134  -4.244  1.00 14.00           C
ATOM      4  O   SER A  42       9.426   7.215  -3.746  1.00 15.75           O
ATOM      5  OG  SER A  42      12.866   5.322  -5.039  1.00 18.50           O
";

#[test]
fn atoms_keep_their_names_and_say_what_they_are() {
    let s = read(SER, Format::Pdb);
    assert_eq!(s.atoms[1].name(), "CA");
    assert_eq!(s.atoms[4].name(), "OG");
    let card = s.card(4, false).unwrap();
    assert_eq!((card.title.as_str(), card.tag.as_str()), ("Serine", "SER 42"));
    assert_eq!(card.chain.as_deref(), Some("A"));
    assert_eq!(card.atom.as_deref(), Some("OG"));
    assert_eq!(card.element.map(|e| e.0), Some("O"));
    assert_eq!(card.b, Some(("B-factor", 18.5)));
    assert!(card.pos.is_some());
    // A residue picked off the ribbon is not one atom.
    let card = s.card(1, true).unwrap();
    assert_eq!(card.ss, Some(crate::protein::Ss::Coil));
    assert!(card.atom.is_none() && card.pos.is_none() && card.element.is_none());
}

#[test]
fn a_lone_molecule_names_its_atoms_by_number() {
    let s = read(WATER_XYZ, Format::Xyz);
    let card = s.card(0, false).unwrap();
    assert_eq!(card.tag, "atom 1");
    assert_eq!(card.element.map(|e| e.0), Some("O"));
}

#[test]
fn long_names_are_cut_at_four_bytes() {
    assert_eq!(&crate::protein::atom_name("HD21X"), b"HD21");
    assert_eq!(&crate::protein::atom_name(" C1'"), b"C1'\0");
}

#[test]
fn a_haem_pocket_holds_its_histidines() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples/4hhb.pdb.gz");
    let s = load_bytes("4hhb.pdb.gz", std::fs::read(path).unwrap()).unwrap();
    let haems = s.copies("HEM");
    assert_eq!(haems.len(), 4);
    let site = s.site(haems[0], 5.0);
    let names: Vec<String> = site.iter().map(|&r| format!("{}{}", s.residues[r].name, s.residues[r].seq)).collect();
    // The alpha chain's proximal and distal histidines hold the haem.
    assert!(names.contains(&"HIS87".to_string()), "{names:?}");
    assert!(names.contains(&"HIS58".to_string()), "{names:?}");
    assert!(site.iter().all(|&r| !s.residues[r].outside_polymer()));
    assert!(s.side_chain_centre(site[0]).is_some());
}

/// The first three nucleotides of the Dickerson dodecamer, 1BNA.
const CGC: &str = "ATOM      1  O5'  DC A   1      18.935  34.195  25.617  1.00 64.35           O  \nATOM      2  C5'  DC A   1      19.130  33.921  24.219  1.00 44.69           C  \nATOM      3  C4'  DC A   1      19.961  32.668  24.100  1.00 31.28           C  \nATOM      4  O4'  DC A   1      19.360  31.583  24.852  1.00 37.45           O  \nATOM      5  C3'  DC A   1      20.172  32.122  22.694  1.00 46.72           C  \nATOM      6  O3'  DC A   1      21.350  31.325  22.681  1.00 48.89           O  \nATOM      7  C2'  DC A   1      18.948  31.223  22.647  1.00 30.88           C  \nATOM      8  C1'  DC A   1      19.231  30.482  23.944  1.00 36.58           C  \nATOM      9  N1   DC A   1      18.070  29.661  24.380  1.00 40.51           N  \nATOM     10  C2   DC A   1      18.224  28.454  25.015  1.00 16.62           C  \nATOM     11  O2   DC A   1      19.360  28.014  25.214  1.00 27.75           O  \nATOM     12  N3   DC A   1      17.143  27.761  25.377  1.00 20.55           N  \nATOM     13  C4   DC A   1      15.917  28.226  25.120  1.00 34.72           C  \nATOM     14  N4   DC A   1      14.828  27.477  25.444  1.00 40.31           N  \nATOM     15  C5   DC A   1      15.719  29.442  24.471  1.00 30.78           C  \nATOM     16  C6   DC A   1      16.843  30.171  24.101  1.00 25.90           C  \nATOM     17  P    DG A   2      22.409  31.286  21.483  1.00 58.85           P  \nATOM     18  OP1  DG A   2      23.536  32.157  21.851  1.00 57.82           O  \nATOM     19  OP2  DG A   2      21.822  31.459  20.139  1.00 78.33           O  \nATOM     20  O5'  DG A   2      22.840  29.751  21.498  1.00 40.36           O  \nATOM     21  C5'  DG A   2      23.543  29.175  22.594  1.00 47.19           C  \nATOM     22  C4'  DG A   2      23.494  27.709  22.279  1.00 47.81           C  \nATOM     23  O4'  DG A   2      22.193  27.252  22.674  1.00 38.76           O  \nATOM     24  C3'  DG A   2      23.693  27.325  20.807  1.00 28.58           C  \nATOM     25  O3'  DG A   2      24.723  26.320  20.653  1.00 40.44           O  \nATOM     26  C2'  DG A   2      22.273  26.885  20.416  1.00 21.14           C  \nATOM     27  C1'  DG A   2      21.721  26.304  21.716  1.00 33.95           C  \nATOM     28  N9   DG A   2      20.237  26.470  21.780  1.00 34.00           N  \nATOM     29  C8   DG A   2      19.526  27.584  21.429  1.00 36.47           C  \nATOM     30  N7   DG A   2      18.207  27.455  21.636  1.00 32.37           N  \nATOM     31  C5   DG A   2      18.083  26.212  22.142  1.00 15.06           C  \nATOM     32  C6   DG A   2      16.904  25.525  22.545  1.00 11.88           C  \nATOM     33  O6   DG A   2      15.739  25.916  22.518  1.00 21.30           O  \nATOM     34  N1   DG A   2      17.197  24.279  23.037  1.00 15.44           N  \nATOM     35  C2   DG A   2      18.434  23.717  23.155  1.00  9.63           C  \nATOM     36  N2   DG A   2      18.508  22.456  23.668  1.00 16.69           N  \nATOM     37  N3   DG A   2      19.537  24.360  22.770  1.00 30.98           N  \nATOM     38  C4   DG A   2      19.290  25.594  22.274  1.00 18.56           C  \nATOM     39  P    DC A   3      25.064  25.621  19.252  1.00 44.67           P  \nATOM     40  OP1  DC A   3      26.506  25.316  19.220  1.00 53.89           O  \nATOM     41  OP2  DC A   3      24.559  26.412  18.115  1.00 57.79           O  \nATOM     42  O5'  DC A   3      24.260  24.246  19.327  1.00 35.42           O  \nATOM     43  C5'  DC A   3      24.584  23.285  20.335  1.00 45.75           C  \nATOM     44  C4'  DC A   3      23.523  22.233  20.245  1.00 43.02           C  \nATOM     45  O4'  DC A   3      22.256  22.844  20.453  1.00 36.85           O  \nATOM     46  C3'  DC A   3      23.424  21.557  18.903  1.00 40.14           C  \nATOM     47  O3'  DC A   3      24.121  20.309  18.928  1.00 49.62           O  \nATOM     48  C2'  DC A   3      21.930  21.406  18.661  1.00 53.79           C  \nATOM     49  C1'  DC A   3      21.278  21.966  19.909  1.00 22.18           C  \nATOM     50  N1   DC A   3      20.196  22.889  19.521  1.00 25.44           N  \nATOM     51  C2   DC A   3      18.909  22.584  19.816  1.00 19.81           C  \nATOM     52  O2   DC A   3      18.685  21.512  20.382  1.00 29.92           O  \nATOM     53  N3   DC A   3      17.935  23.447  19.502  1.00 21.59           N  \nATOM     54  C4   DC A   3      18.217  24.603  18.897  1.00 14.01           C  \nATOM     55  N4   DC A   3      17.221  25.499  18.629  1.00 26.88           N  \nATOM     56  C5   DC A   3      19.526  24.945  18.571  1.00 27.59           C  \nATOM     57  C6   DC A   3      20.537  24.048  18.899  1.00 27.05           C  \n";

#[test]
fn a_dna_strand_is_traced_and_its_bases_are_no_ligands() {
    let s = read(CGC, Format::Pdb);
    assert_eq!(s.residues.len(), 3);
    assert!(s.has_polymer(), "a strand alone is a polymer");
    // The first has no phosphorus, and is traced through its C4'.
    let traced: Vec<&str> = s.residues.iter().map(|r| s.atoms[r.nucleic.unwrap() as usize].name()).collect();
    assert_eq!(traced, ["C4'", "P", "P"]);
    assert!(s.ligands().is_empty());
    assert!(s.residues.iter().all(|r| !r.outside_polymer()));
}

#[test]
fn a_free_nucleotide_stays_a_ligand() {
    // One residue on its own, linked to nothing.
    let one: String = CGC.lines().filter(|l| l.get(22..26).is_some_and(|n| n.trim() == "2")).map(|l| format!("{}\n", l.replacen("ATOM  ", "HETATM", 1))).collect();
    let s = read(&one, Format::Pdb);
    assert!(s.residues[0].nucleic.is_none());
}

#[test]
fn a_helix_sits_where_helices_do_on_the_ramachandran_plot() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples/1crn.pdb.gz");
    let s = load_bytes("1crn.pdb.gz", std::fs::read(path).unwrap()).unwrap();
    let t = s.torsions();
    // The ends have no neighbour to measure across.
    assert!(t[0].is_none() && t[45].is_none());
    let helix: Vec<(f32, f32)> = (0..s.residues.len()).filter(|&r| s.residues[r].ss == crate::protein::Ss::Helix).filter_map(|r| t[r]).collect();
    assert!(helix.len() > 10);
    let (phi, psi) = helix.iter().fold((0.0, 0.0), |(a, b), (p, q)| (a + p, b + q));
    let (phi, psi) = (phi / helix.len() as f32, psi / helix.len() as f32);
    assert!((-80.0..-45.0).contains(&phi) && (-60.0..-25.0).contains(&psi), "{phi} {psi}");
}

#[test]
fn predicted_aligned_error_is_read_as_alphafold_writes_it() {
    let text = r#"[{"predicted_aligned_error":[[0,1,4],[1,0,2.5],[5,3,0]],"max_predicted_aligned_error":31.75}]"#;
    let pae = crate::protein::Pae::from_json(text).unwrap();
    assert_eq!((pae.n, pae.max), (3, 31.75));
    assert_eq!(pae.at(1, 2), 2.5);
    assert_eq!(pae.at(2, 0), 5.0);
    assert!(crate::protein::Pae::from_json(r#"[{"predicted_aligned_error":[[0,1],[1]]}]"#).is_err());
}

#[test]
fn force_field_names_are_amino_acids() {
    for (name, code) in [("HSD", 'H'), ("HIE", 'H'), ("CYX", 'C'), ("ASH", 'D'), ("GLH", 'E'), ("LYN", 'K')] {
        assert_eq!(crate::protein::three_to_one(name), code, "{name}");
    }
}
