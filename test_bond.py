import pytest

from bond import Bond


def test_bond():
    with pytest.raises(ValueError):
        Bond(0, 0)

    b = Bond(0, 1)
    assert b.get_bond_order() == 1

    b1 = Bond(1, 0)
    assert b == b1
