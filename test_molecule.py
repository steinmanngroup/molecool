import atom
from molecule import Molecule, OBMolecule
from util import LABEL2Z

import openbabel

def load_molecule_from_xyz(filename):
    mol = Molecule()
    with open(filename, 'r') as f:
        n = int(f.readline())
        title = f.readline().strip()
        mol.set_name(title)
        for i in range(n):
            tokens = f.readline().split()
            Z = LABEL2Z[tokens[0]]
            mol.add_atom(atom.Atom(Z, xyz=list(map(float, tokens[1:])), idx=i))

    return mol

def test_basic_molecule():
    mol = Molecule()
    assert mol.get_num_atoms() == 0

    # testing basic stuff
    mol.set_name("test")
    assert mol.get_name() == "test"

    mol.set_charge(-5)
    assert mol.get_charge() == -5

    mol.set_multiplicity(1)
    assert mol.get_multiplicity() == 1


def test_molecule_from_file():
    # test with a molecule
    mol = load_molecule_from_xyz('HOH.xyz')
    assert mol.get_num_atoms() == 3
    assert len(list(mol.get_bonds())) == 2
    assert mol.get_name() == "HOH"


def test_molecule_copy():
    mol = load_molecule_from_xyz('CH3.xyz')
    mol2 = Molecule.from_molecule(mol)
    assert mol.get_num_atoms() == mol2.get_num_atoms()
    assert mol.get_name() == mol2.get_name()
    assert mol.get_charge() == mol2.get_charge()
    assert mol.get_multiplicity() == mol2.get_multiplicity()

    # get atoms
    assert mol.get_atom(0) == mol2.get_atom(0)
    assert len(list(mol.get_atoms())) == len(list(mol2.get_atoms()))
    for a1, a2 in zip(mol.get_atoms(), mol2.get_atoms()):
        assert a1 == a2

    # get bonds and angles
    assert len(list(mol.get_bonds())) == len(list(mol2.get_bonds()))
    assert len(list(mol.get_angles())) == len(list(mol2.get_angles()))

    # find_children
    ch1 = mol.find_children(mol.get_atom(0))
    ch2 = mol2.find_children(mol2.get_atom(0))
    assert len(ch1) == len(ch2)
    for a1, a2 in zip(ch1, ch2):
        assert a1 == a2

def test_molecule_copy_openbabel():
    """ tests if copy to openbabel molecule gives same result as regular molecule """
    mol = load_molecule_from_xyz('CH3.xyz')
    mol2 = OBMolecule.from_molecule(mol)
    assert mol.get_num_atoms() == mol2.get_num_atoms()
    assert mol.get_name() == mol2.get_name()
    assert mol.get_charge() == mol2.get_charge()
    assert mol.get_multiplicity() == mol2.get_multiplicity()

    # get atoms
    assert mol.get_atom(0) == mol2.get_atom(0)
    assert len(list(mol.get_atoms())) == len(list(mol2.get_atoms()))
    for a1, a2 in zip(mol.get_atoms(), mol2.get_atoms()):
        assert a1 == a2

    # get bonds and angles
    assert len(list(mol.get_bonds())) == len(list(mol2.get_bonds()))
    assert len(list(mol.get_angles())) == len(list(mol2.get_angles()))

    # find_children
    ch1 = mol.find_children(mol.get_atom(0))
    ch2 = mol2.find_children(mol2.get_atom(0))
    assert len(ch1) == len(ch2)
    for a1, a2 in zip(ch1, ch2):
        assert a1 == a2

def test_openbabel_to_molecule_copy():
    _obmol = OBMoleculeFromFilenameAndFormat('HOH.xyz', file_format='xyz')
    mol = load_molecule_from_xyz('HOH.xyz')
    mol2 = OBMolecule(_obmol)

    assert mol.get_num_atoms() == mol2.get_num_atoms()
    assert mol.get_name() == mol2.get_name()
    assert mol.get_charge() == mol2.get_charge()
    assert mol.get_multiplicity() == mol2.get_multiplicity()

    # get atoms
    assert mol.get_atom(0) == mol2.get_atom(0)
    assert len(list(mol.get_atoms())) == len(list(mol2.get_atoms()))
    for a1, a2 in zip(mol.get_atoms(), mol2.get_atoms()):
        assert a1 == a2

    # get bonds and angles
    assert len(list(mol.get_bonds())) == len(list(mol2.get_bonds()))
    assert len(list(mol.get_angles())) == len(list(mol2.get_angles()))

    # find_children
    ch1 = mol.find_children(mol.get_atom(0))
    ch2 = mol2.find_children(mol2.get_atom(0))
    assert len(ch1) == len(ch2)
    for a1, a2 in zip(ch1, ch2):
        assert a1 == a2

def OBMoleculeFromFilenameAndFormat(filename, file_format='pdb'):
    """ Loads a molecule into an OpenBabel molecule.

        Arguments:
        filename -- the file to load
        file_format -- file format to load

        Returns:
        OpenBabel OBMol instance of the molecule
    """
    obc = openbabel.OBConversion()
    obc.SetInFormat(file_format)
    mol = openbabel.OBMol()
    obc.ReadFile(mol, filename)
    return mol
