class Angle(object):
    """ An angle between three Atom objects.

        Currently, there is no reference to the actual Atom objects but instead
        they are indexed in the parent molecule class.

        An angle between three atoms:

        A
        |
        |
        B----C

        is defined in the Angle class as

        >>> a = Angle(B.getIdx(), A.getIdx(), C.getIdx())

        where B is the vertex of the angle.
    """
    def __init__(self, vertex, id1, id2):
        self._vertex = vertex
        self._id1 = id1
        self._id2 = id2

    def __eq__(self, other):
        this = set([self._vertex, self._id1, self._id2])
        that = set([other._vertex, other._id1, other._id2])
        return len(this - that) == 0

    def __repr__(self):
        return("Angle({0:d},{1:d},{2:d})".format(self._vertex, self._id1, self._id2))
