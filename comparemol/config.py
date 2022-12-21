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
