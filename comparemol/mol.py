from functools import cached_property

import numpy as np

from numpy.typing import ArrayLike


class Mol:
    """Molecule.
    
    Parameters
    ----------
    types : ArrayLike
        Types of atoms, which can be int, str, or any other types.
    coord : ArrayLike
        Coordinates of atoms.

    Examples
    --------
    >>> mol1 = Mol([0, 0], [[0., 0., 0.], [1.,1.,1.]])
    >>> mol2 = Mol([0, 0], [[1., 1., 1.], [2.73205081, 1., 1.]])
    >>> mol1 == mol2
    True
    """
    def __init__(self, types: ArrayLike, coord: ArrayLike) -> None:
        self.coord = np.asarray(coord).reshape(-1, 3)
        self.types = np.asarray(types).reshape(-1)

    @cached_property
    def distance_matrix(self) -> np.ndarray:
        """Distance matrix of the molecule."""
        return np.linalg.norm(self.coord[:, np.newaxis] - self.coord, axis=-1)

    @cached_property
    def sorted_idx(self) -> np.ndarray:
        """Sorted indexes.
        
        Atom indexes sorted in this way:
        (1) Sort by atom type;
        (2) For the same atom type, sort by sorted distance to other atoms.
        """
        sorted_distance = np.sort(self.distance_matrix)[:,::-1]
        features = [x.ravel() for x in np.split(sorted_distance, sorted_distance.shape[1], axis=1)]
        return np.lexsort((*features, self.types))

    @cached_property
    def sorted_distance_matrix(self) -> np.ndarray:
        """Sorted distance matrix of the molecule."""
        return self.distance_matrix[self.sorted_idx][:, self.sorted_idx]

    def __eq__(self, other: "Mol") -> bool:
        """Check if two molecules are equal."""
        return (np.all(self.types[self.sorted_idx] == other.types[other.sorted_idx])
            and np.allclose(self.sorted_distance_matrix, other.sorted_distance_matrix))
