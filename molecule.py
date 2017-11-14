import copy
import numpy

import atom
import bond
import angle

import pytest

class BaseMolecule(object):
    pass

class Molecule(object):
    """ A molecule

        A molecule is a collection of atoms.

        The molecule class can identify all bonds and angles in the molecule.
        This can be quite costly since we use a brute force approach.
    """
    _bond_threshold = 0.45 # Added threshold for bonds. Replicates openbabel

    def __init__(self):
        """ Initializes an empty molecule """
        self._charge = 0
        self._multiplicity = 1
        self._atoms = []
        self._bonds = []
        self._name = ""

    # class methods
    @classmethod
    def fromMolecule(cls, m):
        """ Loads a molecule from another molecule by copying over all data """
        M = cls()
        M.setCharge(m.getCharge())
        M.setMultiplicity(m.getMultiplicity())
        M.setName(m.getName())
        if m.getNumAtoms() > 0:
            M.addAtoms(*m.getAtoms())

        # currently we do not transfer bond information
        return M

    #
    # getters and setters for various properties
    #
    def addAtom(self, _atom):
        #assert isinstance(_atom, atom.Atom), "You attempted to add something that was not an atom."
        self._atoms.append(copy.deepcopy(_atom))

    def addAtoms(self, *args):
        for _atom in args:
            self.addAtom(_atom)

    def getNumAtoms(self):
        """ Returns the number of atoms in the molecule """
        return len(self._atoms)

    def getAtom(self, idx):
        n_atoms = self.getNumAtoms()
        if idx < 0:
            raise IndexError("argument idx to getAtom must be >= 0")
        if idx >= n_atoms:
            raise IndexError("argument idx to getAtom must be < {}".format(n_atoms))
        return self._atoms[idx]

    def getAtoms(self):
        """ Returns all atoms (as an iterator) in the molecule """
        for _atom in self._atoms:
            yield _atom

    def getBonds(self):
        """ Returns all bonds (as an iterator) in the molecule

            If the bond list has not been calculated before, the bonds are
            percieved through the percieveBonds method
        """
        if len(self._bonds) == 0:
            self._bonds = list(self.percieveBonds())

        for _bond in self._bonds:
            yield _bond

    def getName(self):
        return self._name

    def setName(self, value):
        assert isinstance(value, str)
        self._name = value

    def getCharge(self):
        return self._charge

    def setCharge(self, value):
        assert isinstance(value, int)
        self._charge = value

    def getMultiplicity(self):
        return self._multiplicity

    def setMultiplicity(self, value):
        assert isinstance(value, int)
        self._multiplicity = value

    #
    # properties that are lazily evaluated such as bonds and angles
    #
    def percieveBonds(self):
        """ This method attempts to percieve bonds

            It works by comparing atom distances to covalent radii of the atoms.
            It is not optimized in any way.
        """

        for iat, atom1 in enumerate(self.getAtoms()):
            for jat, atom2 in enumerate(self.getAtoms()):
                if iat <= jat: continue
                dr = atom2.getCoordinate() - atom1.getCoordinate()
                R2 = dr.dot(dr)

                dr_cov = atom1.getCovalentRadius() + atom2.getCovalentRadius() + self._bond_threshold
                R2_cov = dr_cov**2
                if R2 < R2_cov:
                    yield bond.Bond(id1=atom1.getIdx(), id2=atom2.getIdx())


    def percieveAngles(self):
        """ This method attemps to percieve angles

            It works by iterating through all bonds in the molecule
        """
        for ibd, bond1 in enumerate(self.getBonds()):
            for jbd, bond2 in enumerate(self.getBonds()):
                if ibd <= jbd: continue
                jatm = bond1.sharesAtom(bond2)
                if jatm >= 0:
                    iatm = bond1.getNbrAtomIdx(jatm)
                    katm = bond2.getNbrAtomIdx(jatm)
                    yield angle.Angle(iatm, jatm, katm)

    #
    # specialized functions to extract information stored in
    # other classes related to molecule
    #
    def getCoordinates(self):
        """ Returns a numpy array with all the coordinates
            of all the atoms in the molecule in angstrom
        """
        c = numpy.zeros((self.getNumAtoms(), 3))
        for iat, _atom in enumerate(self.getAtoms()):
            c[iat] = _atom.getCoordinate()

        return c

    def setCoordinates(self, c):
        """ Sets the coordinates of all atoms in the molecule from
            the numpy array
        """
        assert isinstance(c, numpy.ndarray)
        (n,k) = numpy.shape(c)
        assert n == self.getNumAtoms()
        for iat, _atom in enumerate(self.getAtoms()):
            _atom.setCoordinate(c[iat])


    def getCenterOfMass(self):
        """ Calculates the center of mass of the molecule """
        mass = 0.0
        Rcm = numpy.zeros(3)
        for atom in self.getAtoms():
            atom_mass = atom.getMass()
            mass += atom_mass
            Rcm += atom.getCoordinate() * atom_mass

        assert mass != 0.0, "Total mass of molecule cannot be zero."

        return Rcm / mass

    #
    # Specialized iterators
    #
    def iterAtomAtoms(self, atom):
        """ Iterates over all atoms covalently bound to an atom

            Arguments:
            atom_index -- the atom whose neighbours to get
        """

        neighbour_indices = []
        for ibond, bond in enumerate(self.getBonds()):
            try:
                nbr_index = bond.getNbrAtomIdx(atom.getIdx())
            except ValueError:
                continue
            else:
                neighbour_indices.append(nbr_index)

        for _atom in self._atoms:
            if _atom.getIdx() in neighbour_indices:
                yield _atom


    def findChildren(self, atom):
        """ Finds all atoms which are connected in the molecular graph to atom

            Arguments:
            atom -- the atom to use as the base for finding neighbours
        """
        atoms = []
        old_atoms = [atom]

        for i in range(2):
            atoms = []
            for _atom in old_atoms:
                new_atoms = []
                for _nbr in self.iterAtomAtoms(_atom):
                    if _nbr not in old_atoms:
                        new_atoms.append(_nbr)
                ll = [a.getIdx() for a in new_atoms]

                if len(new_atoms) > 0:
                    atoms.extend(new_atoms)

            if len(atoms) > 0:
                old_atoms.extend(atoms)
            else:
                break

        return old_atoms

def load_molecule(filename):
    from util import LABEL2Z
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

def test_molecule():
    mol = Molecule()
    assert mol.getNumAtoms() == 0

    # testing basic stuff
    mol.setName("test")
    assert mol.getName() == "test"

    mol.setCharge(-5)
    assert mol.getCharge() == -5

    mol.setMultiplicity(1)
    assert mol.getMultiplicity() == 1

    # test with a molecule
    mol = load_molecule('HOH.xyz')
    assert mol.getNumAtoms() == 3
    assert len(list(mol.getBonds())) == 2
    assert mol.getName() == "HOH"

if __name__ == '__main__':
    pass
