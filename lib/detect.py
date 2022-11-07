#!/usr/bin/env python

import sys
import pickle
import numpy as np
from datetime import datetime as dt, timedelta as td

from . import dynfor
from . import docutil
from . import thermodyn

# Take over the contents of dynfor.detect to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.detect, 'detect', sys.modules[__name__])



def cyclone_by_contour(msl, grid):
    ''' Cyclone masks by finding the outermost closed contour
   
    Reimplementation of the Wernli and Schwierz (2006) algorithm, including the
    modifications described in Sprenger et al. (2017). To avoid the technically 
    difficult contour tracing, we base the detection on sorted sea-level 
    pressure values. 
    
    The criteria used in the original implementation of the algorithm are
     - Location: Deepest sea-level pressure minimum within 750 km radius
     - Topography below 1500 m at location
     - Cyclone mask: Enclosed in contour of min 100 km and max 7500 km length
     - No sea-level pressure maximum in the enclosed contour
    
    We adapt the minimum and maximum contour length to cyclone area thresholds
    assuming circular contours. The thresholds then correspond to a minimum 
    size of 800 km^2 and a maximum size of 4.5e6 km^2.
    
    Looping from minimum to maximum values, cyclone masks are then grown by deciding 
    whether each new grid point (1) belongs to an existing cyclone, (2) separates to 
    existing cyclones, (3) constitutes a new sea-level pressure minimum, or (4) none
    of the above.
    
    These cases can be separated using the following rules:
     (1) The point is adjacent to exactly one existing cyclone mask, which does not
         exceed its maximum size
     (2) The point is adjacent to two or more existing cyclone masks
     (3) The point is adjacent to no existing cyclone masks, but constitutes at a local
         minimum
     (4) All points not fulfilling any of the above
    
    In addition to the contour tracking, we further avoid a predefined contour
    interval at which the analysis is carried out. The chosen contour interval however
    implies a minimum prominence of the sea-level pressure minima, which is here
    explicity defined and enforced. Following Wernli and Schwierz (2006) contour
    interval, we set this minimum prominence to 2 hPa.
   
    The algorithm uses the following configration parameters
     - `dynfor.config.cyc_minsize`: Minimum size in km^2
     - `dynfor.config.cyc_maxsize`: Maximum size in km^2
     - `dynfor.config.cyc_maxoro`: Maximum orographic height in the units of
       the `oro` parameter to this function.
     - `dynfor.config.cyc_mindist`: Minimum distance between two cyclone centres in km
     - `dynfor.config.cyc_minprominence`: Minimum prominence of the sea-level
       pressure minimum in the unit of the `msl`/`msls` parameters to this function.
   
    Parameters
    ----------
   
    msl : np.ndarray with shape (nz,ny,nx) and dtype float64
        Sea-level pressure (or any other suitable field).
    grid : gridlib.grid
        Grid and static information about the array. Make sure grid.oro contains the 
        model orography in masl.
   
    Returns
    -------
    np.ndarray with shape (nz,ny,nx) and dtype int32
        Cyclone mask with different integers designating different cyclones
    np.ndarray with shape (nz,nn,5)
        Meta data about each cyclone: (1) latitutde and (2) longitude of its centre,
        (3) minimum SLP, (4) SLP at the outermost contour, and (5) cyclone size.
    '''
    
    # Preprocessing
    if len(msl.shape) == 4:
        if msl.shape[1] > 1:
            raise ValueError('Only one vertical level can be used for the the cyclone detection.')
        msl = msl[:,0,:,:]

    nt, ny, nx = msl.shape
    
    # Creating array holding only the grid index for the respective locations for the following sort
    xidx = np.empty(msl.shape[1:], dtype='i4')
    yidx = np.empty(msl.shape[1:], dtype='i4')
    xidx[:,:] = np.arange(nx)[np.newaxis,:]
    yidx[:,:] = np.arange(ny)[:,np.newaxis]
    
    # Allocating the sorted arrays
    sortshape = (msl.shape[0], msl.shape[1]*msl.shape[2])
    xidx_sort = np.empty(sortshape, dtype='i4')
    yidx_sort = np.empty(sortshape, dtype='i4')
    msl_sort = np.empty(sortshape)

    # Sorting the input
    for k in range(nt):
        sortidx = np.unravel_index(np.argsort(msl[k,:,:], axis=None), msl[k,:,:].shape)
        msl_sort[k,:] = msl[k,:,:][sortidx]
        xidx_sort[k,:] = xidx[sortidx]
        yidx_sort[k,:] = yidx[sortidx]
        
    # The actual cyclone detections
    mask, meta = cyclone_by_contour_fortran(1000, msl, msl_sort, xidx_sort, yidx_sort, grid.oro,
            grid.x[0,:], grid.y[:,0], grid.dx, grid.dy)

    return mask, meta


def block_by_grad_rev(dat, grid, lat_band=(30,70),
        local_move=(27,36), total_move=(40,54), min_duration=21):
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
    bi = dynfor.detect.block_indicator_grad_rev(dat, grid.dx, grid.dy)

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
            blockmask[tidx,j0:j1,:] = dynfor.utils.mask_minimum_connect(
                    bi[tidx,:,:], 
                    np.array(seeds[tidx], dtype='i4'), 
                    0.0
            )

    blockmask = blockmask.astype('i1')

    return blockmask


def frontalvolume_largescale(tfp, dx, dy, mountain_mask=None):
    ''' Detect frontal zones as coherent areas with strong TFP gradients
    
    The large-scale version of this function (applicable for example to ERA-Interim data)
    uses only the large-scale gradient threshold. 

    The gradient threshold can be set by ``dynfor.config.tfp_mingrad_largescale``,
    the minimum size of the frontal zones by ``dynfor.config.tfp_grad_minsize``.

    Parameters
    ----------

    tfp : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
        Thermal front parameter (TFP) field, typically potential or equivalent potential 
        temperature.
    dx : np.ndarray with shape (ny,nx) and dtype float64
        The double grid spacing in x-direction to be directly for centered differences.
        ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``. 
    dy : np.ndarray with shape (ny,nx) and dtype float64
        The double grid spacing in y-direction to be directly for centered differences.
        ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
    mountain_mask : np.ndarray with dtype bool and shape (nt,nz,ny,nx) or (ny,nx)
        *Optional*, default ``None``. The gradient of tfp will be set to 0 at masked grid points.
    
    Returns
    -------
    np.ndarray with shape (nz,ny,nx) and dtype int32
        Detected frontal zones labeled with object IDs > 0, zero everywhere else.

    See also
    --------
    :meth:`frontalvolume_smallscale`, :meth:`frontline`
    '''

    thres = dynfor.config.tfp_mingrad_largescale        # 4.5e-5
    min_size = dynfor.config.minsize            # 75000.0e+6
    cellsize = np.abs(dx)*np.abs(dy)/4.0

    labels = np.empty(tfp.shape, dtype='i4')
    for tidx in range(tfp.shape[0]):
        sys.stdout.write('%d/%d\r' % (tidx+1, tfp.shape[0]))
        sys.stdout.flush()
        ddx, ddy = dynfor.derivatives.grad(tfp[tidx,:,:,:], dx, dy)
        tfp_grad = np.sqrt(ddx**2 + ddy**2)
        if not type(mountain_mask) == type(None):
            if len(mountain_mask.shape) == 4:
                tfp_grad[mountain_mask[tidx,:,:,:]] = 0.
            elif len(mountain_mask.shape) == 2:
                tfp_grad[:,mountain_mask] = 0.
            else:
                raise ValueError('Mountain mask must be either 2-dimensional (y,x) or 4-dimensional (t,z,y,x).')
            
        mask = tfp_grad > thres

        labels_, sizes = dynfor.utils.label_connected_3d(mask, cellsize, 2500)
        for n, size in zip(range(1,len(sizes)+1), sizes):
            if size == 0:
                break

            if size < (mask.shape[0]*min_size):
                labels_[labels_ == n] = 0

        labels[tidx,:,:,:] = labels_
    
    print('')

    return labels


def frontalvolume_smallscale(tfp, dx, dy, mountain_mask=None, tfps=None, maxobj=3000, quiet=True):
    ''' Detect frontal zones as coherent areas with strong TFP gradients
    
    The small-scale version of this function (applicable for example to NORA10 data)
    uses a combination of a large-scale gradient threshold for a smoothed TFP field
    and a small-scale gradient threshold for the unfiltered TFP field.

    The gradient thresholds can be set by ``dynfor.config.tfp_mingrad_largescale``,
    and ``dynfor.config.tfp_mingrad_smallscale``, the minimum size of the frontal 
    zones by ``dynfor.config.tfp_grad_minsize`` and the amount of smoothing by
    ``dynfor.config.nsmooth``.

    Parameters
    ----------

    tfp : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
        Thermal front parameter (TFP) field, typically potential or equivalent potential 
        temperature.
    dx : np.ndarray with shape (ny,nx) and dtype float64
        The double grid spacing in x-direction to be directly for centered differences.
        ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``. 
    dy : np.ndarray with shape (ny,nx) and dtype float64
        The double grid spacing in y-direction to be directly for centered differences.
        ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
    mountain_mask : np.ndarray with shape (nt,nz,ny,nx) and dtype bool
        *Optional*, default ``None``. The gradient of tfp will be set to 0 at masked grid points 
    tfps: np.ndarray with shape (nt,nz,ny,nx) and dtype float64
        Smoothed thermal front parameter (TFP) field. Must be created out of the same tfp field
        as is passed as first argument
        *optional*, default ``None``. If ``None``, tfps is created by dynfor.utils.smooth_xy
    maxobj : int
        *Option*, default 3000. Maximum number of objects to be found.
    quiet : bool
        *Optional*, default ``True``. If ``False`` display progress by current time step.
    
    Returns
    -------
    np.ndarray with shape (nz,ny,nx) and dtype int32
        Detected frontal zones labeled with object IDs > 0, zero everywhere else.

    See also
    --------
    :meth:`frontalvolume_largescale`, :meth:`frontline`
    '''

    nsmooth = dynfor.config.nsmooth
    thres_ls = dynfor.config.tfp_mingrad_largescale        # 4.5e-5
    thres_ss = dynfor.config.tfp_mingrad_smallscale        # 7.5e-5
    min_size = dynfor.config.minsize            # 75000.0e+6
    cellsize = np.abs(dx)*np.abs(dy)/4.0

    labels = np.empty(tfp.shape, dtype='i4')
    for tidx in range(tfp.shape[0]):
        if not quiet:
            print(f'{tidx}/{tfp.shape[0]}', end=chr(13))
        ddx, ddy = dynfor.derivatives.grad(tfp[tidx,:,:,:], dx, dy)
        tfp_grad = np.sqrt(ddx**2 + ddy**2)
        if not type(mountain_mask) == type(None):
            tfp_grad[mountain_mask[tidx,:,:,:]] = 0.
            
        if type(tfps) == type(None):
            tfps_ = dynfor.utils.smooth_xy(tfp[tidx,:,:,:], nsmooth)
        else:
            tfps_ = tfps[tidx,:,:,:]
            
        ddx, ddy = dynfor.derivatives.grad(tfps_, dx, dy)
        tfps_grad = np.sqrt(ddx**2 + ddy**2)
        if not type(mountain_mask) == type(None):
            tfps_grad[mountain_mask[tidx,:,:,:]] = 0.

        mask = (tfp_grad > thres_ss) & (tfps_grad > thres_ls)

        labels_, sizes = dynfor.utils.label_connected_3d(mask, cellsize, maxobj)
        for n, size in zip(range(1,len(sizes)+1), sizes):
            if size == 0:
                break

            if size < (mask.shape[0]*min_size):
                labels_[labels_ == n] = 0

        labels[tidx,:,:,:] = labels_

    if not quiet:
        print('Finished.')

    return labels


def atmospheric_river(mflux_u, mflux_v, grid, thresh_ivt=250,thresh_dist=2000, merthresh=None):
    '''
    Atmospheric river detection after Rutz et. al (2014)
    DOI: 10.1175/MWR-D-13-00168.1

    Atmospheric rivers are defined as continous areas of at least 2000 km length with 
    an atmospheric water vapour transport of 250 kg m-1 s-1 (both these threshold can 
    be changed)

    Instead of the standard 250 kg m-1 s-1 threshold a variable threshold can be used (e.g. the 85 percentile)
    Moreover, compared to Rutz et al. (2014) a meridional transport threshold of 50 kg m-1 s-1 is used.

    Parameters
    ----------
    mflux_u : np.ndarray with shape (ny,nx)
        vertical integral of eastward water vapour flux
    mflux_v : np.ndarray with shape (ny,nx)
        vertical integral of northward water vapour flux
    grid :  gridlib.grid
        Grid information as provided by metopen or get_instantaneous.
    thresh_ivt : float
        threshold used in detecting high ivt areas if another variable is used (like tcw) this threshold should be changed. 
        Default : 250 kg m-1 s-1
    thresh_dist : float
        Minimum length (in km) for the area required to be defined as an atmospheric river. Default: 2000
    merthresh : float 
        Meridional transport threshold. Default : False (but common threshold is 50) 
    
    Returns
    -------
    riv_mask : 
        Mask with 1 at the grid points where the atmospheric rivers are detected, and 0 otherwise
    riv_nr : 
        Mask in which each AR is associated with a different number. 
        
    '''
    from scipy.signal import convolve2d
    from scipy.ndimage import label,generate_binary_structure
    from dynlib.utils import dist_sphere
    #Define output arrays
    riv_mask = np.zeros(mflux_u.shape) 
    riv_nr   = np.zeros(mflux_v.shape)

    #Thresh_ivt should be either one value or an array of similar size as the ivt array

    #1a. Calculate areas with ivt > threshold
    ivt = (mflux_u**2 + mflux_v**2)**0.5
    ivt_mask = (ivt > thresh_ivt)*1.0

    #1b. Nr. the respective areas
    ivt_areas, num_features = label(ivt_mask) #Nr. Areas

    #Correct for areas across the West/Eastern boundaries
    for y in range(ivt_areas.shape[0]):
        if ivt_areas[y, 0] > 0 and ivt_areas[y, -1] > 0:
            ivt_areas[ivt_areas == ivt_areas[y, -1]] = ivt_areas[y, 0]

    #2. Trace the contour by calculating the nr. of neighbours with ivt > threshold
    ivt_mask_conv = convolve2d(ivt_mask,np.ones([3,3]),mode="same",boundary="wrap")

    #3. Loop over different areas to check if maximum distance > distance threshold
    rivnrtemp = 1
    for idx in range(1,num_features+1): 

        templats = grid.y[(ivt_areas == idx) & (ivt_mask_conv < 8.0)]
        templons = grid.x[(ivt_areas == idx) & (ivt_mask_conv < 8.0)]

        #In step 1b some of the areas have been merged, resulting in that
        #Sometimes an area is not existing anymore
        if(len(templats) == 0):
            continue 

        #Boolean to assign if an area is an atmospheric river or not
        riv = False

        #Check if object straddles the equator
        maxlat = np.nanmax(templats)
        minlat = np.nanmin(templats)
        if((maxlat > 0) & (minlat < 0)):
            continue

        #Check if mean poleward transport is above threshold
        if not type(merthresh) ==type(None): 
            mean_mflux_v = np.nanmean(mflux_v[(ivt_areas == idx)])
            if(np.abs(mean_mflux_v) < merthresh):
                continue

        #Calculate maximum difference in lat and lon (in degrees)
        difflat = np.abs(maxlat - minlat) #np.abs(np.nanmin(templats) - np.nanmax(templats))
        difflon = np.abs(np.nanmin(templons) - np.nanmax(templons))
        if(difflon > 180.0):
            difflon = np.abs(difflon -360.0)
        #First check on maximum distance (Check if maximum possible distance < treshold)
        if(((difflat**2 + difflon**2)**0.5)*158 < thresh_dist):
            continue
        #First check if latitude difference is > threshold (then it is an atmospheric river)
        elif(difflat*111.0 > thresh_dist):
            #print("Difflat > thresh")
            riv = True
        #Else: check great circle distance fo all outermost points of area
        else:
            #print("Checking with great circle distance")
            for idx1 in range(len(templats)):
                for idx2 in range(len(templats)-1,-1,-1):  
                    #Calculate great circle distance
                    if(idx1 != idx2):
                        dist = dist_sphere(templats[idx1], templons[idx1], templats[idx2], templons[idx2])
                        #If distance is > threshold distance, stop checking
                        if(dist > thresh_dist):
                            riv= True 
                            break
                if(riv == True):
                    break
        #If the area is an atmospheric river, save to output
        if(riv == True):
            riv_mask[(ivt_areas == idx)] = 1.0
            riv_nr[(ivt_areas == idx)]   = rivnrtemp
            rivnrtemp += 1

    #Return: mask with areas belonging to an atmospheric river (riv_mask)
    #And for every river area a specific nr. 
    return riv_mask, riv_nr


def cold_air_outbreak_index(t850, msl, sst, ci, lsm, cao_thres=None):
    ''' Calculate a cold air outbreak-index, optionally apply threshold
    
    The index is defined by the difference in potential temperature between the sea-surface
    and 850 hPa. The function returns this difference field where larger than zero if no
    threshold is given (default), otherwise the function returns a mask field marking regions
    where the threshold is exceeded.

    Land-regions are masked by NaN. Grid points are regarded as land if the land-sea mask 
    exceeds 0.5.
    
    Parameters
    ----------

    t850 : np.ndarray with shape (nt,ny,nx) and dtype float64
        Temperature field at 850 hPa.
    msl : np.ndarray with shape (nt,ny,nx) and dtype float64
        Sea-level pressure field.
    sst : np.ndarray with shape (nt,ny,nx) and dtype float64
        Sea-surface temperature field.
    ci : np.ndarray with shape (nt,ny,nx) and dtype float64
        Sea-ice concentration field.
    lsm : np.ndarray with shape (ny,nx) and dtype float64
        Land-sea mask, invariant in time, where 1 marks land.
    
    Returns
    -------
    np.ndarray with shape (nz,ny,nx) and dtype float64 or bool
        Cold air outbreak index or detected cold air outbreak events.
    '''

    kappa = dynfor.consts.rl / dynfor.consts.cp
    landmask = lsm > 0.5

    # Convert to theta_850 and theta_sfc
    t850 *= (1000/850)**kappa
    tsfc = thermodyn.theta_from_temp(sst, msl)
    caoidx = tsfc - t850 

    # Mask land and sea-ice
    caoidx[...,landmask] = np.nan
    caoidx[ci > 0.5] = np.nan
    
    # Apply threshold or mask only negative CAO index values
    if not type(cao_thres) == type(None):
        caoidx = caoidx > cao_thres
    else:
        caoidx[caoidx < 0] = 0.0
    
    return caoidx



#
