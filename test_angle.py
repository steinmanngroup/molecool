import pytest

from angle import Angle
def test_angle():
    with pytest.raises(ValueError):
        Angle(0, 0, 1)

    with pytest.raises(ValueError):
        Angle(-1, 1, 2)

    with pytest.raises(ValueError):
        Angle(1, -1, 2)

    with pytest.raises(ValueError):
        Angle(1, 1, -2)

    a1 = Angle(0, 1, 2)
    a2 = Angle(0, 2, 1)
    assert a1 == a2

    a3 = Angle(1, 2, 0)
    assert a3 != a1
    assert a3 != a2

    a4 = Angle(0, 1, 3)
    assert a4 != a1
