import numpy as np

from comparemol import Mol, get_rotatation


def test_apply():
    mol1 = Mol([0, 0], [[1., 1., 1.], [2., 2., 2.]])
    mol2 = Mol([0, 0], [[0., 0., 0.], [1.73205081, 0., 0.]])
    r = get_rotatation(mol1, mol2)
    assert np.allclose(r.apply([1.73205081, 0., 0.]), np.ones(3))
    assert mol1 == mol2
