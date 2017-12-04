import util

class Formatter(object):
    pass

class XYZFormatter(Formatter):
    def __init__(self):
        Formatter.__init__(self)

    @classmethod
    def coordinate(cls, atom):
        return "{0:<2s}{1[0]:20.9f}{1[1]:16.9f}{1[2]:16.9f}\n".format(util.Z2LABEL[atom.get_nuclear_charge()], atom.get_coordinate())

    @classmethod
    def coordinates(cls, atoms):
        s = ""
        for atom in atoms:
            s += cls.coordinate(atom)
        return s
