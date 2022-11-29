from comparemol import Mol

def test_translate():
    mol1 = Mol([0, 0], [[0., 0., 0.], [1., 1., 1.]])
    mol2 = Mol([0, 0], [[1., 1., 1.], [2., 2., 2.]])
    assert mol1 == mol2


def test_rotate():
    mol1 = Mol([0, 0], [[0., 0., 0.], [1., 1., 1.]])
    mol2 = Mol([0, 0], [[0., 0., 0.], [1.73205081, 0., 0.]])
    assert mol1 == mol2


def test_exchange():
    mol1 = Mol([0, 0], [[0., 0., 0.], [1., 1., 1.]])
    mol2 = Mol([0, 0], [[1., 1., 1.], [0., 0., 0.]])
    assert mol1 == mol2


def test_exchange_atype():
    mol1 = Mol([0, 0, 1], [[0., 0., 0.], [1., 1., 1.], [2., 2., 2.]])
    mol2 = Mol([0, 1, 0], [[0., 0., 0.], [1., 1., 1.], [2., 2., 2.]])
    print(mol1.sorted_idx)
    print(mol2.sorted_idx)
    assert mol1 != mol2


def test_not_equal():
    mol1 = Mol([0, 0], [[0., 0., 0.], [1., 1., 1.]])
    mol2 = Mol([0, 0], [[1., 1., 1.], [1., 2., 2.]])
    assert mol1 != mol2
