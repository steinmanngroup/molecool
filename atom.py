import numpy

import util

class Atom(object):
    """ An atom

        The minimum amount of information required is the
        nuclear charge Z of the atom.

        A lot of additional default information is added to the atom
        based on the nuclear charge if not provided by the instantia-
        tor. The information is:

          * mass
          * Van der Waal radius
          * covalent radius
          * coordination number
          * label

        Arguments:
        Z -- nuclear charge of atom. This argument is mandatory.

        Keyword Arguments:
        mass -- the mass of the atom in atomic units. Default is specified by using the nuclear charge.
        xyz -- the Cartesian coordinate of the atom in Angstrom. Default is origo.
        idx -- the atom index. If not specified, a value of -1 is assigned.
        fcharge -- formal charge used for cat- and anions. Default 0.
    """
    def __init__(self, Z, **kwargs):
        assert Z > 0, "Nuclear charge of atom must be greater than zero."
        self._z = Z
        self._c = numpy.array(kwargs.get('xyz', [0, 0, 0]))
        self._mass = kwargs.get('mass', util.MASSES[Z])
        self._idx = kwargs.get('idx', -1)
        self._fcharge = kwargs.get('fcharge', 0)
        self._vdw_radius = kwargs.get('vwdradius', util.VDWRADII[Z])
        self._cov_radius = kwargs.get('covradius', util.COVALENTRADII[Z])
        self._coordination = kwargs.get('coordination', util.COORDINATION[Z])
        self._label = util.Z2LABEL[Z]

    @classmethod
    def fromAtom(cls, other):
        atom = cls(other.getNuclearCharge(),
                   mass=other.getMass(),
                   fcharge=other.getFormalCharge(),
                   xyz=other.getCoordinate(),
                   idx=other.getIdx())
        return atom

    @classmethod
    def fromOBAtom(cls, _obatom):
        """ Constructs an Atom from an openbabel OBAtom """
        x, y, z = _obatom.GetX(), _obatom.GetY(), _obatom.GetZ()
        _atom = cls(_obatom.GetAtomicNum(),
                    xyz=[x, y, z],
                    idx=_obatom.GetId(), # not using internal index (GetIdx) from openbabel because it gets remapped
                    fcharge=_obatom.GetFormalCharge())
        return _atom

    def getMass(self):
        """ Returns the mass of the atom """
        return self._mass

    def getNuclearCharge(self):
        """ Returns the nuclear charge of the atom """
        return self._z

    def getFormalCharge(self):
        """ Returns the formal charge of the atom """
        return self._fcharge

    def setFormalCharge(self, value):
        """ Sets the integer formal charge of the atom

            Arguments:
            value -- the integer formal charge
        """
        if not isinstance(int, value):
            raise TypeError

        self._fcharge = value

    def getIdx(self):
        """ Returns the internal index of the atom """
        return self._idx

    def setIdx(self, value):
        """ Sets the internal index of the atom

            Arguments:
            value -- the integer index
        """
        if not isinstance(int, value):
            raise TypeError
        self._idx = value

    def getLabel(self):
        """ Returns the human readable atomic label """
        return self._label

    def getCoordinate(self):
        """ returns the 3D coordinate of the atom """
        return self._c

    def setCoordinate(self, value):
        """ Sets the 3D coordinate of the atom

            Arguments:
            value -- the coordinate given as a numpy array
        """
        if not isinstance(value, numpy.ndarray):
            raise TypeError("Argument 'value' must be of type numpy array")
        (n, ) = numpy.shape(value)
        if n != 3:
            raise ValueError("Dimensions of data do not match. Expected 3 but got {}".format(n))
        self._c = value

    def getVDWRadius(self):
        """ Returns the Van der Waal radius of the atom """
        return self._vdw_radius

    def getCovalentRadius(self):
        """ Returns the covalent radius of the atom """
        return self._cov_radius

    def setCoordination(self, value):
        """ Sets the coordination number

            Arguments:
            value -- the coordination number
        """
        if not isinstance(int, value):
            raise TypeError

        max_coordination = util.COORDINATION[self._z]
        if value > max_coordination:
            raise ValueError("Coordination number too large.")

        self._coordination = value

    def getCoordination(self):
        """ Returns the coordination number """
        return self._coordination

    def __eq__(self, other):
        """ Tests that two atoms are equal using various measures """
        EPS=1.0e-6
        dr = self.getCoordinate() - other.getCoordinate()
        R2 = dr.dot(dr)

        return self.getIdx() == other.getIdx() and numpy.sqrt(R2) < EPS

    def __repr__(self):
        return "Atom({0:d}, xyz=[{1[0]:.7f}, {1[1]:.7f}, {1[2]:.7f}, idx={2:d}])".format(self.getNuclearCharge(), self.getCoordinate(), self.getIdx())
