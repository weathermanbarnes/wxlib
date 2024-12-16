##### Import packages

import sys
import pickle
import numpy as np
import scipy.interpolate as interp
from datetime import datetime as dt, timedelta as td

##### Import created packages
from utils import utils

#####
def block_by_grad_rev(dat, grid, lat_band=(30,70),
        local_move=(27,36), total_move=(40,54), min_duration=21,
        block_dj=10, block_di=10, grid_cyclic_ew=True,threshold=0.0,**kwargs):

    from blocking_hepers import block_indicator_grad_rev

    ''' Detect blocking following the Masato et al. (2012) procedure based on persistent gradient reversals

    The procedure involves three main steps:
     
     1. Find instantaneous regions with reversed gradients (done by 
        :meth:`block_indicator_grad_rev`)
     2. Connect these regions in time and apply stationarity criteria
     3. Mask those regions with reversed gradients that fullfil the criteria in step 2
    
    Parameters
    ----------

    dat : np.ndarray with dimensions (nt,ny,nx)
        Data to base the detection on, e.g. 500 hPa geopotential or PV on 330 K. Blocks
        should be associated with a positive anomaly of the respective field, so for 
        example Northern Hemisphere PV requires a sign switch.
    grid : dynlib.gridlib.grid
        Grid information as provided by metopen or get_instantaneous.
    lat_band: 2-tuple of float or int
        Which latitudes in degrees to consider for the blocking detection. Default:
        ``(30, 70)``. Both latitudes must coincide with a grid latitude.
    local_move: 2-tuple of float or int
        Maximum movement in degrees latitude/longitude of the block center from 
        one time step to the next. Default: ``(27, 36)``.
    total_move: 2-tuple of float or int
        Maximum movement in degrees latitude/longitude of the block center during 
        the entire duration. Default: ``(40, 54)``.
    min_duration: int
        Minimum number of time steps for the block to persist. Default: ``21``.

    Returns
    -------
    np.ndarray of dtype i1 and dimensions (nt,ny,nx)
        Binary mask field indicating detected blocks.
    '''

    import scipy.ndimage.filters as ndf


    print('Stage 1+2: Finding gradient reversals and conenct them in time')
    ny, nx = grid.x.shape
    tlen = dat.shape[0]

    # Stage 1: Calculate blocking indicator
    bi = block_indicator_grad_rev(dat, grid.dx, grid.dy, grid.nx, grid.ny, grid.nt)

    # Translate thresholds from degrees to grid point indexes
    j0 = np.argwhere(grid.y[:,0] == lat_band[0])[0,0]
    j1 = np.argwhere(grid.y[:,0] == lat_band[1])[0,0]
    if j0 > j1: 
        j1, j0 = j0, j1

    dx = abs(sorted(grid.x[0,1:]-grid.x[0,:-1])[nx//2])
    dy = abs(sorted(grid.y[1:,0]-grid.y[:-1,0])[ny//2])

    djl_max, dil_max = local_move[0]/dy, local_move[1]/dx
    djt_max, dit_max = total_move[0]/dy, total_move[1]/dx

    # Stage 2: Connect in time
    bi = bi[:,j0:j1,:]
    
    # Identify local maxima
    extrema = (ndf.maximum_filter(bi, size=(1,3,3), mode='wrap') == bi)
    extrema[bi <= 0] = False

    # Remove local extrema directly at the border of the considered domain
    extrema[:,0,:] = False
    extrema[:,-1,:] = False

    # Joining into block objects
    blocks = []
    prevblocks = {}
    for tidx, date in enumerate(grid.t_parsed):
        done = {}
        curblocks = {}
        # From mask to grid point indexes
        expos = np.argwhere(extrema[tidx,:,:] > 0)
        for prevpos, blocknr in prevblocks.items():
            initpos = blocks[blocknr]['pos'][0]
            dij2 = 9999999
            for pos in expos:
                pos = tuple(pos)
                # Movement since last
                djl = abs(pos[0]-prevpos[0])
                dil = abs(pos[1]-prevpos[1])
                if dil > nx/2:
                    dil = nx - dil
                # Total movement
                djt = abs(pos[0]-initpos[0])
                dit = abs(pos[1]-initpos[1])
                if dit > nx/2:
                    dit = nx - dit
                # Check if connected to known block
                if djl <= djl_max and djt <= djt_max and \
                   dil <= dil_max and dit <= dit_max:
                    # New closest continuation of the block
                    if djl**2 + dil**2 < dij2:
                        dij2 = djl**2 + dil**2
                        savepos = pos
                    # Belongs to same block structure, but another maximum is closer
                    else:
                        pass

                    done[pos] = True

            if dij2 < 9999999:
                curblocks[savepos] = blocknr
                blocks[blocknr]['pos'].append(savepos)
                blocks[blocknr]['blockidx'].append(bi[tidx,savepos[0],savepos[1]])

        # Save new blocks
        for pos in expos:
            pos = tuple(pos)
            if not pos in done:
                blocks.append({
                    'pos': [pos, ],
                    'blockidx': [bi[tidx,pos[0],pos[1]], ],
                    'onset': date,
                })
                curblocks[pos] = len(blocks) -1

        prevblocks = curblocks

    #return blocks,bi
    
    #def blocks_part2():
    # Create an list of seeds.
    print('Stage 3: Applying minimum duration criterion and write out block masks')
    seeds = []
    for tidx in range(tlen):
        seeds.append([])
    for block in blocks:
        if len(block['blockidx']) >= min_duration:
            tidx0 = list(grid.t_parsed).index(block['onset'])
            for dtidx in range(len(block['blockidx'])):
                pos = block['pos'][dtidx]
                # convert to Fortran indexes
                seeds[tidx0+dtidx].append([pos[0]+1, pos[1]+1])
    
    blockmask = np.zeros(dat.shape, dtype='bool')
    for tidx in range(bi.shape[0]):
        if len(seeds[tidx]) > 0:
            blockmask[tidx,j0:j1,:] = utils.find_regions_above_threshold(
                    bi[tidx,:,:], 
                    np.array(seeds[tidx], dtype='i4'), 
                    threshold=10e-6
            )

    blockmask = blockmask.astype('i1')

    return blockmask