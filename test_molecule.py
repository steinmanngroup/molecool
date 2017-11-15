import atom
from molecule import Molecule
from util import LABEL2Z

def load_molecule_from_xyz(filename):
    mol = Molecule()
    with open(filename, 'r') as f:
        n = int(f.readline())
        title = f.readline().strip()
        mol.setName(title)
        for i in range(n):
            tokens = f.readline().split()
            Z = LABEL2Z[tokens[0]]
            mol.addAtom(atom.Atom(Z, xyz=map(float, tokens[1:]), idx=i))

    return mol

def test_basic_molecule():
    mol = Molecule()
    assert mol.getNumAtoms() == 0

    # testing basic stuff
    mol.setName("test")
    assert mol.getName() == "test"

    mol.setCharge(-5)
    assert mol.getCharge() == -5

    mol.setMultiplicity(1)
    assert mol.getMultiplicity() == 1


def test_molecule_from_file():
    # test with a molecule
    mol = load_molecule_from_xyz('HOH.xyz')
    assert mol.getNumAtoms() == 3
    assert len(list(mol.getBonds())) == 2
    assert mol.getName() == "HOH"


def test_molecule_copy():
    mol = load_molecule_from_xyz('CH3.xyz')
    mol2 = Molecule.fromMolecule(mol)
    assert mol.getNumAtoms() == mol2.getNumAtoms()
    assert mol.getName() == mol2.getName()
    assert mol.getCharge() == mol2.getCharge()
    assert mol.getMultiplicity() == mol2.getMultiplicity()

    # get atoms
    assert mol.getAtom(0) == mol2.getAtom(0)
    assert len(list(mol.getAtoms())) == len(list(mol2.getAtoms()))
    for a1, a2 in zip(mol.getAtoms(), mol2.getAtoms()):
        assert a1 == a2

    # get bonds and angles
    assert len(list(mol.getBonds())) == len(list(mol2.getBonds()))
    assert len(list(mol.getAngles())) == len(list(mol2.getAngles()))

    # findChildren
    ch1 = mol.findChildren(mol.getAtom(0))
    ch2 = mol2.findChildren(mol2.getAtom(0))
    assert len(ch1) == len(ch2)
    for a1, a2 in zip(ch1, ch2):
        assert a1 == a2
