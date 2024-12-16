##### Import packages

import sys
import pickle
import numpy as np
import scipy.interpolate as interp
from datetime import datetime as dt, timedelta as td

##### Import created packages
from utils import utils

#####

def block_indicator_grad_rev(dat, dx, dy, nx, ny, nz, 
                             block_dj=10, block_di=10, grid_cyclic_ew=True,
                             **kwargs):
    """
    Calculates the local block intensity indicator field using large-scale gradient reversals.

    Parameters
    ----------
    dat : np.ndarray of shape (nz, ny, nx)
        PV (or any other suitable field) used to detect gradient reversals.
    dx : np.ndarray of shape (ny, nx)
        The double grid spacing in x-direction for centered differences.
        dx[j, i] is the x-distance between (j, i+1) and (j, i-1).
    dy : np.ndarray of shape (ny, nx)
        The double grid spacing in y-direction for centered differences.
        dy[j, i] is the y-distance between (j+1, i) and (j-1, i).
    nx : int
        Grid size in x-direction.
    ny : int
        Grid size in y-direction.
    nz : int
        Grid size in z- or t-direction.
    block_dj : int
        The half-width for local blocking index.
    block_di : int
        The half-width for longitudinal smoothing.
    grid_cyclic_ew : bool
        If True, the grid is cyclic in the east-west direction.

    Returns
    -------
    res : np.ndarray of shape (nz, ny, nx)
        Local block intensity indicator field. Maxima in this field are blocking centers 
        that can be tracked in time and space to verify the minimum persistence and 
        stationarity requirements.
    """

    # Initialize output array
    res = np.zeros((nz, ny, nx), dtype=np.float64)
    tmp = np.zeros((nz, ny, nx), dtype=np.float64)

    for k in range(nz):
        # Local blocking index
        for j in range(ny - block_dj, block_dj, -1):
            for i in range(nx):
                tmp[k, j, i] = (np.sum(dat[k, j - block_dj:j, i]) - np.sum(dat[k, j + 1:j + block_dj + 1, i])) / block_dj
        
        # Longitudinal smoothing by running mean
        for j in range(ny - block_dj, block_dj, -1):
            for i in range(1 + block_di, nx - block_di):
                res[k, j, i] = np.sum(tmp[k, j, i - block_di:i + block_di + 1]) / (2 * block_di + 1)
    
    # Handle periodic grid for longitudinal smoothing
    if grid_cyclic_ew:
        for k in range(nz):
            for j in range(ny - block_dj, block_dj, -1):
                for i in range(1, block_di + 1):
                    res[k, j, i] = (np.sum(tmp[k, j, nx - (block_di - i):nx]) + 
                                    np.sum(tmp[k, j, 0:i + block_di])) / (2 * block_di + 1)
                
                for i in range(nx - block_di, nx):
                    res[k, j, i] = (np.sum(tmp[k, j, i - block_di:nx]) + 
                                    np.sum(tmp[k, j, 0:i - (nx - block_di)])) / (2 * block_di + 1)

    return res