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
