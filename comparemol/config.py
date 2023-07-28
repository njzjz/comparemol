"""Configurations for tolerances."""
from typing import Tuple

rtol = 1e-4
atol = 1e-4


def set_tol(rt: float, at: float):
    """Set the tolerance for the comparison of two molecules.

    Parameters
    ----------
    rt : float
        Relative tolerance.
    at : float
        Absolute tolerance.

    Examples
    --------
    >>> set_tol(1e-5, 1e-5)
    """
    global rtol, atol
    rtol = rt
    atol = at


def get_tol() -> Tuple[float, float]:
    """Get the tolerance for the comparison of two molecules.

    Returns
    -------
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance.

    Examples
    --------
    >>> get_tol()
    (1e-05, 1e-05)
    """
    return rtol, atol
