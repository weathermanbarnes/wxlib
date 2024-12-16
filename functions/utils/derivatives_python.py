import numpy as np

def grad(dat, dx, dy, grid_cyclic_ew=True):
    """
    Calculates second partial derivative in x and y directions using 2nd-order centered differences.
    It returns NaN on the first and last latitude and longitude for non-cyclic grids.

    Parameters
    ----------
    dat : np.ndarray of shape (nz, ny, nx)
        Data array.
    dx : np.ndarray of shape (ny, nx)
        The double grid spacing in x-direction, directly for centered differences.
        dx[j, i] is the x-distance between (j, i+1) and (j, i-1).
    dy : np.ndarray of shape (ny, nx)
        The double grid spacing in y-direction, directly for centered differences.
        dy[j, i] is the y-distance between (j+1, i) and (j-1, i).
    grid_cyclic_ew : bool
        If True, the grid is cyclic in the east-west direction.

    Returns
    -------
    resx : np.ndarray of shape (nz, ny, nx)
        x-derivative of the input data.
    resy : np.ndarray of shape (nz, ny, nx)
        y-derivative of the input data.
    """

    # Get dimensions
    nz, ny, nx = dat.shape
    
    # Initialize output arrays
    resx = np.full((nz, ny, nx), np.nan)
    resy = np.full((nz, ny, nx), np.nan)
    
    # Calculate x-derivative
    resx[:, 1:ny-1, 1:nx-1] = (dat[:, 1:ny-1, 2:nx] - dat[:, 1:ny-1, 0:nx-2]) / dx[1:ny-1, 1:nx-1]
    
    # Calculate y-derivative
    resy[:, 1:ny-1, 1:nx-1] = (dat[:, 2:ny, 1:nx-1] - dat[:, 0:ny-2, 1:nx-1]) / dy[1:ny-1, 1:nx-1]
    
    # Set to NaN where gradient results are not valid: j = 1 and ny
    resx[:, 0, :] = np.nan
    resy[:, 0, :] = np.nan
    resx[:, ny-1, :] = np.nan
    resy[:, ny-1, :] = np.nan
    
    if not grid_cyclic_ew:
        # Set to NaN where gradient results are not valid: i = 1 and nx
        resx[:, :, 0] = np.nan
        resy[:, :, 0] = np.nan
        resx[:, :, nx-1] = np.nan
        resy[:, :, nx-1] = np.nan

    return resx, resy

def lap2(dat, dx, dy, grid_cyclic_ew=True):
    """
    Calculates 2D Laplacian of dat using 2nd-order centered differences.
    Returns NaN on the first and last latitude and longitude for non-cyclic grids.

    Parameters
    ----------
    dat : np.ndarray of shape (nz, ny, nx)
        Data array.
    dx : np.ndarray of shape (ny, nx)
        The double grid spacing in x-direction, directly for centered differences.
        dx[j, i] is the x-distance between (j, i+1) and (j, i-1).
    dy : np.ndarray of shape (ny, nx)
        The double grid spacing in y-direction, directly for centered differences.
        dy[j, i] is the y-distance between (j+1, i) and (j-1, i).
    grid_cyclic_ew : bool
        If True, the grid is cyclic in the east-west direction.

    Returns
    -------
    res : np.ndarray of shape (nz, ny, nx)
        Laplacian of the input data.
    """

    # Get dimensions
    nz, ny, nx = dat.shape

    # Initialize output array
    res = np.full((nz, ny, nx), np.nan)

    # Main calculation for non-boundary points
    res[:, 1:ny-1, 1:nx-1] = (
        (dat[:, 1:ny-1, 2:nx] + dat[:, 1:ny-1, 0:nx-2] - 2 * dat[:, 1:ny-1, 1:nx-1]) / dx[1:ny-1, 1:nx-1]**2 +
        (dat[:, 2:ny, 1:nx-1] + dat[:, 0:ny-2, 1:nx-1] - 2 * dat[:, 1:ny-1, 1:nx-1]) / dy[1:ny-1, 1:nx-1]**2
    )

    if grid_cyclic_ew:
        # Handle cyclic conditions in x-direction (east-west)
        res[:, 1:ny-1, 0] = (
            (dat[:, 1:ny-1, 1] + dat[:, 1:ny-1, nx-1] - 2 * dat[:, 1:ny-1, 0]) / dx[1:ny-1, 0]**2 +
            (dat[:, 2:ny, 0] + dat[:, 0:ny-2, 0] - 2 * dat[:, 1:ny-1, 0]) / dy[1:ny-1, 0]**2
        )
        res[:, 1:ny-1, nx-1] = (
            (dat[:, 1:ny-1, 0] + dat[:, 1:ny-1, nx-2] - 2 * dat[:, 1:ny-1, nx-1]) / dx[1:ny-1, nx-1]**2 +
            (dat[:, 2:ny, nx-1] + dat[:, 0:ny-2, nx-1] - 2 * dat[:, 1:ny-1, nx-1]) / dy[1:ny-1, nx-1]**2
        )
    else:
        # Set boundaries in x-direction to NaN
        res[:, 1:ny-1, 0] = np.nan
        res[:, 1:ny-1, nx-1] = np.nan

    # Set boundaries in y-direction to NaN
    res[:, 0, :] = np.nan
    res[:, ny-1, :] = np.nan

    return res

