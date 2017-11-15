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
        return self._mass

    def getNuclearCharge(self):
        return self._z

    def getFormalCharge(self):
        return self._fcharge

    def setFormalCharge(self, value):
        if not isinstance(int, value):
            raise TypeError

        self._fcharge = value

    def getIdx(self):
        return self._idx

    def setIdx(self, idx):
        if not isinstance(int, idx):
            raise TypeError
        self._idx = value

    def getLabel(self):
        return self._label

    def getCoordinate(self):
        return self._c

    def setCoordinate(self, value):
        (n, ) = numpy.shape(value)
        assert n == 3, "Dimensions of data do not match. Expected 3 but got {}".format(n)
        self._c = value

    def getVDWRadius(self):
        return self._vdw_radius

    def getCovalentRadius(self):
        return self._cov_radius

    def setCoordination(self, value):
        if not isinstance(int, value):
            raise TypeError

        max_coordination = util.COORDINATION[self._z]
        if value > max_coordination:
            raise ValueError("Coordination number too large.")

        self._coordination = value

    def getCoordination(self):
        return self._coordination

    def __eq__(self, other):
        return self.getIdx() == other.getIdx()

    def __repr__(self):
        return "Atom({0:d}, xyz=[{1[0]:.7f}, {1[1]:.7f}, {1[2]:.7f}, idx={2:d}])".format(self.getNuclearCharge(), self.getCoordinate(), self.getIdx())
