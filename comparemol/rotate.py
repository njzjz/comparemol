"""Rotate vector."""
import numpy as np
from scipy.spatial.transform import Rotation as R

from .config import get_tol
from .mol import Mol


def get_rotation(mol1: Mol, mol2: Mol) -> R:
    """Get rotation: mol2 -> mol1.

    Parameters
    ----------
    mol1 : Mol
        molecule 1
    mol2 : Mol
        molecule 2

    Returns
    -------
    scipy.spatial.transform.Rotation
        rotation
    """
    assert mol1 == mol2
    # align the atoms
    coord1 = mol1.coord[mol1.sorted_idx]
    coord2 = mol2.coord[mol2.sorted_idx]
    # align the first atom
    coord1 -= coord1[0]
    coord2 -= coord2[0]
    r, rmse = R.align_vectors(coord1, coord2)
    rtol, atol = get_tol()
    if not np.isclose(rmse, 0.0, rtol=rtol * 10, atol=atol * 10):
        raise RuntimeError(
            f"Molecules not aligned. RMSE={rmse}, rtol={rtol}, atol={atol}"
        )
    return r
