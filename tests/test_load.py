"""Test load from external packages."""
import dpdata
import numpy as np

from comparemol import Mol


def test_load_from_dpdata():
    """Test load from dpdata."""
    mol1 = Mol(["H", "H"], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    system = dpdata.System(
        data={
            "orig": np.zeros(3),
            "atom_names": ["H"],
            "atom_numbs": [2],
            "atom_types": np.zeros(2, dtype=int),
            "coords": np.array([[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]]),
            "cells": np.zeros((1, 3, 3)),
        }
    )
    mol2 = Mol.load_from_dpdata(system)
    assert mol1 == mol2
