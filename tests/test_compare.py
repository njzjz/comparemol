"""Test compare module."""
from comparemol import Mol


def test_translate():
    """Test translate."""
    mol1 = Mol([0, 0], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    mol2 = Mol([0, 0], [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    assert mol1 == mol2


def test_rotate():
    """Test rotate."""
    mol1 = Mol([0, 0], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    mol2 = Mol([0, 0], [[0.0, 0.0, 0.0], [1.73205081, 0.0, 0.0]])
    assert mol1 == mol2


def test_exchange():
    """Test exchange."""
    mol1 = Mol([0, 0], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    mol2 = Mol([0, 0], [[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]])
    assert mol1 == mol2


def test_exchange_atype():
    """Test exchange atom type."""
    mol1 = Mol([0, 0, 1], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    mol2 = Mol([0, 1, 0], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    print(mol1.sorted_idx)
    print(mol2.sorted_idx)
    assert mol1 != mol2


def test_not_equal():
    """Test not equal."""
    mol1 = Mol([0, 0], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    mol2 = Mol([0, 0], [[1.0, 1.0, 1.0], [1.0, 2.0, 2.0]])
    assert mol1 != mol2
