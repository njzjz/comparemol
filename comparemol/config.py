rtol = 1e-4
etol = 1e-4


def set_tol(rt: float, et: float):
    """Set the tolerance for the comparison of two molecules.
    
    Parameters
    ----------
    rt : float
        Relative tolerance.
    et : float
        Absolute tolerance.
    
    Examples
    --------
    >>> set_tol(1e-5, 1e-5)
    """
    global rtol, etol
    rtol = rt
    etol = et
