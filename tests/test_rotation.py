"""Test rotation."""
import numpy as np

from comparemol import Mol, get_rotation


def test_apply():
    """Test apply rotation."""
    mol1 = Mol([0, 0], [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    mol2 = Mol([0, 0], [[0.0, 0.0, 0.0], [1.73205081, 0.0, 0.0]])
    r = get_rotation(mol1, mol2)
    assert np.allclose(r.apply([1.73205081, 0.0, 0.0]), np.ones(3))
    assert mol1 == mol2
