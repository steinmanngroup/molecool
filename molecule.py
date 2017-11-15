import copy
import numpy

import atom
import bond
import angle

import openbabel

class BaseMolecule(object):
    """ A molecule

        A molecule is a collection of atoms and provides methods to
        extract or manipulate atoms.

        Internally, all atomic coordinates are stored in Angstrom
        and this convention should be adhered to when attempting to
        update the coordinates.

        Molecules from other APIs derive from this base class and
        implements all methods defined by this abstract class.
        See for example both Molecule and OBMolecule classes below.

    """
    _bond_threshold = 0.45 # Added threshold for bonds. Replicates openbabel

    def __init__(self):
        """ Initializes an empty molecule """
        self._charge = 0
        self._multiplicity = 1
        self._atoms = []
        self._bonds = []
        self._name = ""

    #
    # class methods
    #
    @classmethod
    def fromMolecule(cls, m):
        """ Loads a molecule from another molecule by copying over all data """
        M = cls()
        M.setCharge(m.getCharge())
        M.setMultiplicity(m.getMultiplicity())
        M.setName(m.getName())
        if m.getNumAtoms() > 0:
            M.addAtoms(*m.getAtoms())

        return M

    #
    # methods to add information to the molecule
    #
    def addAtom(self, _atom):
        """ Adds an atom to the molecule """
        raise NotImplementedError


    def addAtoms(self, *args):
        """ Adds multiple atoms to a molecule """
        for _atom in args:
            self.addAtom(_atom)

    #
    # getters and setters for various properties
    #
    def getNumAtoms(self):
        """ Returns the number of atoms in the molecule """
        raise NotImplementedError


    def getAtom(self, idx):
        """ Returns an atom based on its index

            Arguments:
            idx -- the index (from 0 to getNumAtoms() -1)

            Returns:
            an atom with the specified index
        """
        raise NotImplementedError


    def getAtoms(self):
        """ Returns an iterator of all atoms in the molecule """
        raise NotImplementedError


    def getBonds(self):
        """ Returns an iterator of all bonds in the molecule """
        raise NotImplementedError


    def getAngles(self):
        """ Returns an iterator of all angles in the molecule """
        raise NotImplementedError

    #
    # getters and setters for simple properties
    #
    def getName(self):
        """ Returns the name of the molecule """
        return self._name


    def setName(self, value):
        """ Sets the name of the molecule

            Arguments:
            value -- the name of the molecule
        """
        assert isinstance(value, str)
        self._name = value


    def getCharge(self):
        """ Returns the integer charge of the molecule """
        return self._charge


    def setCharge(self, value):
        """ Sets the integer charge of the molecule

            Arguments:
            value -- the integer charge of the molecule
        """
        if not isinstance(value, int):
            raise ValueError("Argument 'value' must be of type integer")
        self._charge = value


    def getMultiplicity(self):
        """ Returns the multiplicity of the molecule """
        return self._multiplicity


    def setMultiplicity(self, value):
        """ Sets the multiplicity of the molecule

            Arguments:
            value -- the integer charge of the molecule
        """
        if not isinstance(value, int):
            raise ValueError("Argument 'value' must be of type integer")
        self._multiplicity = value

    #
    # properties that are lazily evaluated such as bonds
    #
    def percieveBonds(self):
        """ Attempts to percieve covalent bonds in the molecule """
        raise NotImplementedError

    #
    # specialized getters and setters to extract information stored
    # in other classes related to the molecule such as atoms
    #
    def getCoordinates(self):
        """ Returns a numpy array with all the coordinates of the molecule

            Note: coordinates are always in Angstrom.
        """
        c = numpy.zeros((self.getNumAtoms(), 3))
        for iat, _atom in enumerate(self.getAtoms()):
            c[iat] = _atom.getCoordinate()

        return c


    def setCoordinates(self, value):
        """ Sets coordinates of molecule

            Arguments:
            value -- coordinates to store. Must be numpy array.
        """
        raise NotImplementedError


    def getCenterOfMass(self):
        """ Calculates the center of mass in units of Angstrom """
        mass = 0.0
        Rcm = numpy.zeros(3)
        for _atom in self.getAtoms():
            atom_mass = _atom.getMass()
            mass += atom_mass
            Rcm += _atom.getCoordinate() * atom_mass

        assert mass != 0.0, "Total mass of molecule cannot be zero."

        return Rcm / mass

    #
    # Specialized iterators
    #
    def iterAtomAtoms(self, atom):
        """ Iterates over all atoms covalently bound to an atom

            Arguments:
            atom -- the atom whose neighbours to get
        """

        neighbour_indices = []
        for ibond, bond in enumerate(self.getBonds()):
            try:
                nbr_index = bond.getNbrAtomIdx(atom.getIdx())
            except ValueError:
                continue
            else:
                neighbour_indices.append(nbr_index)

        for _atom in self.getAtoms():
            if _atom.getIdx() in neighbour_indices:
                yield _atom


    def findChildren(self, atom):
        """ Finds all atoms in the molecular graph as the supplied atom

            Arguments:
            atom -- the atom whose neighbours to get

            Returns:
            List of atoms sorted according to their internal index
        """
        raise NotImplementedError


class Molecule(BaseMolecule):
    """ A molecule

        A molecule is a collection of atoms and provides methods to
        extract or manipulate atoms.

        A library-free implementation of the Molecule class. The
        goal is provide a solution (albeit it might be slow) that
        does not depend on any third-party libraries.

        Internally, all atomic coordinates are stored in Angstrom
        and this convention should be adhered to when attempting to
        update the coordinates.
    """
    def __init__(self):
        """ Initializes an empty molecule """
        BaseMolecule.__init__(self)


    def addAtom(self, _atom):
        """ Adds an atom to the molecule """
        self._atoms.append(copy.deepcopy(_atom))


    def getNumAtoms(self):
        """ Returns the number of atoms in the molecule """
        return len(self._atoms)


    def getAtom(self, idx):
        """ Returns an atom based on its index

            Arguments:
            idx -- the index (from 0 to getNumAtoms() -1)

            Returns:
            an atom with the specified index
        """
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
        """ Returns an iterator of all angles in the molecule

            If the bond list has not been calculated before, the bonds are
            percieved through the percieveBonds method
        """
        if len(self._bonds) == 0:
            self._bonds = list(self.percieveBonds())

        for _bond in self._bonds:
            yield _bond

    def getAngles(self):
        """ Returns an iterator of all angles in the molecule """
        for ibd, bond1 in enumerate(self.getBonds()):
            for jbd, bond2 in enumerate(self.getBonds()):
                if ibd <= jbd: continue
                jatm = bond1.sharesAtom(bond2)
                if jatm >= 0:
                    iatm = bond1.getNbrAtomIdx(jatm)
                    katm = bond2.getNbrAtomIdx(jatm)
                    yield angle.Angle(jatm, iatm, katm)


    def percieveBonds(self):
        """ Attempts to percieve covalent bonds in the molecule

            It compares atom distances to covalent radii of the atoms.
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


    def setCoordinates(self, value):
        """ Sets coordinates of molecule

            Arguments:
            value -- coordinates to store. Must be numpy array.
        """
        if not isinstance(value, numpy.ndarray):
            raise TypeError("Argument 'value' must be of type numpy array")
        (n,k) = numpy.shape(value)
        if n != self.getNumAtoms():
            raise ValueError("Argument 'value' has the wrong number of atoms")
        for iat, _atom in enumerate(self.getAtoms()):
            _atom.setCoordinate(value[iat])


    def findChildren(self, atom):
        """ Finds all atoms in the molecular graph as the supplied atom

            Arguments:
            atom -- the atom whose neighbours to get

            Returns:
            List of atoms sorted according to their internal index
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

        return sorted(old_atoms, key=lambda _atom: _atom.getIdx())


class OBMolecule(BaseMolecule):
    """ A molecule

        A molecule is a collection of atoms and provides methods to
        extract or manipulate atoms.

        This is an implementation of the Molecule class that directly
        supports openbabel as the underlying engine.

        Internally, all atomic coordinates are stored in Angstrom
        and this convention should be adhered to when attempting to
        update the coordinates.
    """
    def __init__(self):
        """ Initializes an empty molecule """
        BaseMolecule.__init__(self)
        self._obmol = openbabel.OBMol()


    def addAtom(self, _atom):
        """ Adds an atom to the molecule

            This method creates OBAtoms based on the atoms so it can use
            the internal functions that openbabel provides through OBMol.

            Note: This will surely break at some point for .pdb files etc.
        """
        _obatom = openbabel.OBAtom()
        _obatom.SetAtomicNum(_atom.getNuclearCharge())
        x, y, z = _atom.getCoordinate()
        _obatom.SetVector(x, y, z)
        _obatom.SetId(_atom.getIdx())
        _obatom.SetFormalCharge(_atom.getFormalCharge())
        self._obmol.AddAtom(_obatom)


    def getMultiplicity(self):
        """ Returns the multiplicity of the molecule """
        return self._obmol.GetTotalSpinMultiplicity()


    def getCharge(self):
        """ Returns the integer charge of the molecule """
        return self._obmol.GetTotalCharge()


    def getNumAtoms(self):
        """ Returns the number of atoms in the molecule """
        return self._obmol.NumAtoms()


    def getAtom(self, idx):
        """ Returns an atom based on its index

            Note: The internal format for indices in openbabel is
                  to count from 1 to getNumAtoms() which is offset
                  by one compared to the format we have chosen in
                  the molecool library

            Arguments:
            idx -- the index (from 0 to getNumAtoms() -1)

            Returns:
            an atom with the specified index

        """
        n_atoms = self.getNumAtoms()
        if idx < 0:
            raise IndexError("argument idx to getAtom must be >= 0")
        if idx >= n_atoms:
            raise IndexError("argument idx to getAtom must be < {}".format(n_atoms))
        return atom.Atom.fromOBAtom(self._obmol.GetAtom(idx+1))


    def getAtoms(self):
        """ Returns all atoms (as an iterator) in the molecule """
        for _obatom in openbabel.OBMolAtomIter(self._obmol):
            _atom = atom.Atom.fromOBAtom(_obatom)
            yield _atom


    def getBonds(self):
        """ Returns an iterator of all bonds in the molecule """
        self._obmol.ConnectTheDots() # what happens if called more than once?
        for _obbond in openbabel.OBMolBondIter(self._obmol):
            iat = _obbond.GetBeginAtom()
            jat = _obbond.GetEndAtom()
            id1 = _obbond.GetBeginAtomIdx()
            id2 = _obbond.GetEndAtomIdx()
            yield bond.Bond(iat.GetId(), jat.GetId())


    def getAngles(self):
        """ Returns an iterator of all angles in the molecule """
        self._obmol.FindAngles() # does nothing if angles have been found
        for _obangle in openbabel.OBMolAngleIter(self._obmol):
            vertex, id1, id2 = _obangle
            yield angle.Angle(vertex, id1, id2)


    def percieveBonds(self):
        """ Attempts to percieve covalent bonds in the molecule """
        self._obmol.ConnectTheDots()


    def setCoordinates(self, c):
        """ Sets coordinates of molecule

            Arguments:
            value -- coordinates to store. Must be numpy array.
        """
        if not isinstance(value, numpy.ndarray):
            raise TypeError("Argument 'value' must be of type numpy array")
        (n,k) = numpy.shape(value)
        if n != self.getNumAtoms():
            raise ValueError("Argument 'value' has the wrong number of atoms")
        for iat, _obatom in enumerate(openbabel.OBMolAtomIter(self._obmol)):
            x, y, z = c[iat]
            _obatom.SetVector(x, y, z)


    def findChildren(self, atom):
        """ Finds all atoms in the molecular graph as the supplied atom

            This method uses the internal one supplied with openbabel

            Arguments:
            atom -- the atom whose neighbours to get
        """
        self._obmol.ConnectTheDots() # what happens if called more than once?
        idx = int(atom.getIdx()+1)
        vector = openbabel.vectorInt()
        self._obmol.FindChildren(vector, 0, idx)
        indices = sorted([i-1 for i in vector] + [idx-1])
        atoms = [self.getAtom(i) for i in indices]
        return sorted(atoms, key=lambda _atom: _atom.getIdx())

