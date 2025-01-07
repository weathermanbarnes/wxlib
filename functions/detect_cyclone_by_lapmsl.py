##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
################################ Functions for tracking cyclones by lapmsl (UniMelb Method) ######################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

import sys
import pickle
import numpy as np
import pandas as pd
import scipy.interpolate as interp
from datetime import datetime as dt, timedelta as td
import scipy.interpolate as interp
import scipy.ndimage as ndimg
import cartopy.crs as ccrs

##### Import created packages
from utils import utils, gridlib
from utils import derivatives_python as derivatives

##################################################################################################################################
###################################################### Main functions ############################################################
##################################################################################################################################

def detect_cyclone_by_lapmsl(msl, oro, grid, prev_cyc=None, prev_tracks=None, quiet=False, return_data=True, lmsl_thres_closed=2.0e-09, lmsl_thres_open=7.0e-09, maxdist_lmsl_center=500.0e3, msl_min_prominence_closed=150, msl_min_prominence_open=25, mindist_centers=750.0e3,maxdist_track=750.0e3,tstep=10800,min_lifetime=24,min_dist_travelled=500.0e3, outpath='',outfilename='cyclone_lapmsl_test', **kwargs):
    
    ''' Detect and track cyclone centers following the Melbourne algorithm (MS91)
    DEVELOPED BY: Clemens Spensberger (Bergen)
    Some small adpations have been made to adapt it for wxlib (ARC Coe for 21st Century Weather, Monash University). These include,
    the way in which grid information is given to the software and inclusion of output writing etc. 

    The original detection algorithm is defined in Murray and Simmonds (1991a,b; MS91) developed at UniMelb. This is a simplified
    reimplementation of the detection algorithm complemented by a more fancy tracking. 

    Simplifications compared to MS91:
    - As in MS91, msl is interpolated to polar stereographic grids before the cyclone detection and tracking.
      However, the grid resolution is chosen fine enough that cyclone centeres are simply determined at the
      grid-point resolution.
    - No attempt is made in determining cyclone sizes
    - Prediction of cyclone movement is determined from Kalman-filtered previous cyclone movement rather
      than from the background flow. This requires much fewer assumptions and configuration and reduces the 
      required input data to only sea-level pressure (or equivalent).
    
    This is the (so-far) first function in dynlib to depend on the pandas library. pandas is imported only 
    within this function such that the remainder of dynlib remains usable even when pandas is not available.
    
    Parameters
    ----------
   
    msl : np.ndarray with shape (nz,ny,nx) and dtype float64
        Sea-level pressure (or any other suitable field), pre-filtered to the resolution at which cyclones are
        to be detected.
    grid : gridlib.grid
        Grid and static information about the array. Make sure grid.oro contains the model orography in masl,
        and make sure the orograpghy data is pre-filtered as desired. Typically the same filtering will be 
        applied to grid.oro as to msl.
    prev_cyc : dict of pd.DataFrame
        *Optional*, default ``None``. Cyclone detections for the time step before, separately for the NH and 
        the SH. Required if cyclones detected in this call to the function are meant to continue tracks from 
        a previous call. The return value can be directly passed on towards the next call.
    prev_tracks : pd.Dataframe
        *Optional*, default ``None``. Potentially unfinished cyclone tracks. Required if cyclones detected 
        in this call to the function are meant to continue tracks from a previous call.
    quiet : bool
        *Optional*, default ``False``. If ``True``, progress information is suppressed.
    return_data: bool,
        *Optional*, default ``True``. If ``False``, data is not returned to a written variable. This is useful
        if you only want to write out to a file.
    lmsl_thres_closed : float64
        *Optional*, default ``2.0e-09``. Lower threshold in Pa m**-2 on MSLP laplacian for closed systems.
    lmsl_thres_open : float64
        *Optional*, default ``7.0e-09``. Lower threshold in Pa m**-2 on MSLP laplacian for open systems
    maxdist_lmsl_center : float64
        *Optional*, default ``500.0e3``. Maximum distance in km between maximum in MSLP laplacian and either 
        MSLP minimum (closed system) or grad(MSLP)**2 minimum (open system).
    msl_min_prominence_closed : float64
        *Optional*, default ``150``. Lower threshold in Pa on local prominence compared to environment for a 
        closed system.
    msl_min_prominence_open : float64
        *Optional*, default ``25``. Lower threshold in Pa on local prominence compared to environment for an 
        open system.
    mindist_centers : float64
        *Optional*, default ``500.0e3``. Lower threshold in km on distance between individual identified 
        cyclone centres.
    maxdist_track : float64
        *Optional*, default ``500.0e3``. Upper threshold in km on distance travelled by a cyclone over one 
        time step.
    tstep : float64
        *Optional*, default ``10800``. Time interval in s between two consecutive time steps.
    min_lifetime : float64
        *Optional*, default ``24``. Lower threshold in hours on duration of a tracked cyclone.
    min_dist_travelled : float64
        *Optional*, default ``500.0e3``. Lower threshold in km on total travelled distance during the 
        lifetime of a tracked cyclone.
    outpath: str
        Path to which the files must be written out to. Must end in a "/". Default: ''.
    outfilename: str 
        Filename for both the netcdf and pickle files
   
    Returns (if return_data=True)
    -------
    pd.DataFrame
        Table of all detected cyclone positions regardless of track criteria and associated metadata
    pd.DataFrame
        Table of detected cyclone positions for every point along every track and associated metadata
    pd.DataFrame
        Table of cyclone tracks, including a set of metadata for each track for conventient track filtering
    dict of pd.DataFrame
        Cyclone detections for the last time step to be potentially passed on towards the next call this 
        function, cf. optional parameter ``prev_cyc``.
    pd.DataFrame
        Table of detected cyclone positions including all temporary internal columns belonging to tracks that
        potentially extend beyond the end of the time interval for this call. To be passed on to the next call
        of this function, cf. optional parameter ``prev_tracks``.
    '''
    
    print('Detecting cyclones by lapmsl method ...')
    
    # Prepare grid objects for the polar stereographic grids on which cyclones are being detected and tracked
    timer = dt.now()
    ilons, ilats = utils.concat1lonlat(grid['longitude2D'], grid['latitude2D'])
    ilats = ilats[::-1,0]
    ilons = ilons[0,:]
    ifunc = interp.RectBivariateSpline(ilats, ilons, utils.concat1(oro[::-1,:]))

    if not quiet:
        print(f'   - Sorting grids: {(dt.now()-timer).total_seconds():6.1f} seconds ...')
    maps = {}
    grids = {}
    for hemis in ['nh', 'sh']:
        maps[hemis], grids[hemis] = make_polar_grid(hemis)
        grids[hemis].oro = ifunc(grids[hemis].lat, grids[hemis].lon, grid=False)
    
    timer = dt.now()
    # Initialise loop variables
    if prev_cyc is None:
        cyclones = []
        prev_cyc = {}
        ntracks = 0
    else:
        cyclones = [prev_tracks,]
        ntracks = prev_tracks.track_id.max() + 1
        
    # Loop over time
    for tidx in range(msl.shape[0]):
        if not quiet:
            print(f'   - Detection & tracking: {tidx:03d}/{msl.shape[0]:03d}', end=chr(13))
        ifunc = interp.RectBivariateSpline(ilats, ilons, utils.concat1(msl[tidx,::-1,:]))
    
        for hemis in ['nh', 'sh']:
            hemis_msl = ifunc(grids[hemis].lat, grids[hemis].lon, grid=False)
    
            new_cyclones = locate_cyclones(grid['time'][tidx], hemis_msl, grids[hemis], hemis,\
                                           lmsl_thres_closed, lmsl_thres_open, maxdist_lmsl_center, \
                                            msl_min_prominence_closed, msl_min_prominence_open, \
                                                mindist_centers,)
            cyclones.append(new_cyclones)
    
            if hemis in prev_cyc:
                ntracks = track_cyclones(ntracks, prev_cyc[hemis], new_cyclones, maxdist_track=maxdist_track,tstep=tstep)
    
            prev_cyc[hemis] = new_cyclones
       
    # Concatenate all cyclone detections    
    cyclones = pd.concat(cyclones, ignore_index=True)
    if not quiet:
        print(f'   - Detection & tracking: {(dt.now()-timer).total_seconds():6.1f} seconds.')
    
    timer = dt.now()
    # Now that the basic detection and tracking is complete: run the Kalman smoother for a best estimate
    # of cyclone position and velocity
    tracks = []
    trkids = set(cyclones.track_id)
    
    for trkid in trkids:
        if trkid < 0:
            continue
        track = cyclones[cyclones.track_id == trkid]
        if len(track) > 1:
            smooth_cyclone_track(maps, cyclones, track, tstep)
        track = cyclones[cyclones.track_id == trkid]
        tracks.append(aggregate_trackinfo(track))
    
    tracks = pd.concat(tracks, ignore_index=True)
    if not quiet:
        print(f'   - Smoothing & trackinfo: {(dt.now()-timer).total_seconds():6.1f} seconds.')
    
    timer = dt.now()
    # Filter tracks for minimum length, minimum movement, etc.
    if not quiet:
        print('   - Filtering tracks ...')
    tracks, unfinished_tracks = filter_tracks(tracks, lastdate=grid['time'][-1],min_lifetime=min_lifetime, min_dist_travelled=min_dist_travelled)
    
    sel_track_ids = set(unfinished_tracks.track_id)
    sel_cyclones = np.array([cyc.track_id in sel_track_ids for cycid, cyc in cyclones.iterrows()])
    unfinished_tracks = cyclones.iloc[sel_cyclones].reset_index(drop=True)
    
    sel_track_ids = set(tracks.track_id)
    sel_cyclones = np.array([cyc.track_id in sel_track_ids for cycid, cyc in cyclones.iterrows()])
    tracked_cyclones = cyclones.iloc[sel_cyclones].reset_index(drop=True)
    
    # Clean-up temporary columns for all cyclones belonging to certainly finished tracks
    cyclones = cyclones.drop([
            'Px_pred', 'Py_pred', 'Px_filter', 'Py_filter',
            'x', 'y', 'x_pred', 'y_pred', 'x_filter', 'y_filter', 'x_smooth', 'y_smooth',
            'u', 'v', 'hemis',
        ], axis=1)
    tracked_cyclones = tracked_cyclones.drop([
            'Px_pred', 'Py_pred', 'Px_filter', 'Py_filter',
            'x', 'y', 'x_pred', 'y_pred', 'x_filter', 'y_filter', 'x_smooth', 'y_smooth',
            'u', 'v', 'hemis',
        ], axis=1)
    
    if not quiet:
        print(f'Filtering & cleanup: {(dt.now()-timer).total_seconds():6.1f} seconds.')
    
    if outfilename is not None:
        print('   -  Writing cyclone files ...')
        with open(outpath+outfilename+'_cyclones.pkl', 'wb') as handle:
            pickle.dump(cyclones, handle)
        with open(outpath+outfilename+'_tracked_cyclones.pkl', 'wb') as handle:
            pickle.dump(tracked_cyclones, handle)
        with open(outpath+outfilename+'_tracks.pkl', 'wb') as handle:
            pickle.dump(tracks, handle)
        with open(outpath+outfilename+'_prev_cyc.pkl', 'wb') as handle:
            pickle.dump(prev_cyc, handle)
        with open(outpath+outfilename+'_unfinished_tracks.pkl', 'wb') as handle:
            pickle.dump(unfinished_tracks, handle)
    
    print('Completed!')
    if return_data:
        return cyclones, tracked_cyclones, tracks, prev_cyc, unfinished_tracks

def cyclone_clusters(str_id, str_lon, str_lat, str_dt,
     distthresh=1.0,timthresh=36.0,lngthresh=1.5,timlngthresh=48.0):
    
    '''
    The basis idea of the clustering algorithm is that it checks if multiple cyclone tracks follow a 
    similar path, based on the 'cyclone families' described in Bjerknes and Solberg (1922). For details
    see Weijenborg and Spengler (2024). The algorith further divides cyclone clusters into two different
    types, a 'Bjerknes type' close to the cyclone families of Bjerkens and Solberg (1922) and a stagnant
    type. The former type detects cyclones that follow each other over a certain minimum distance, whereas
    the stagnant type includes cyclones which do not move much in space, but still have a proximity over time.


    Parameters
    ----------
    
    str_id: np.array of length N, with cyclone id number for every point along every track
    str_lat: np.array of length N, with corresponding latitude positions for every point along every track
    str_lon: np.array of length N,  with corresponding longitude positions for every point along every track
    str_dt: np.array of length N, with corresponding time (as .. array) for every point along every track
    
    distthresh = 1.0 #1. Distance criterium (in Rossby Radii)
    timthresh = 36.0 #2. Time criterium (in hours)
    lngthresh = 1.5 #3. Length overlap criterium (in Rossby Radii) 
    timlngthresh = 48.0 #4. Time overlap criterium (in hours)
    
    Returns
    -------
    list
        A list of cyclones which are uniquely clustered together, all cyclones are listed (even 'clusters' 
        of length 1)
    list
        A list of cyclones which are uniquely clustered together according to the Bjerknes type definition
    list
        A list of cyclones which are uniquely clustered together according to the stagnant type definition
    np.array
        np.array of length N, with a 1 if particular point is connected to any other point, for every point 
        along every track
    '''

    from .cluster_helpers import connect_cyclones, find_cluster_type_dokm, get_indices_sparse, unnest
    from scipy.sparse import dok_matrix
    
    #Create options dictionary (easier to pass to different functions)
    Options={'distthresh': distthresh,'timthresh': timthresh,'lngthresh': lngthresh,'timlngthresh': timlngthresh}
    
    #Results array for 
    str_connected   = np.zeros(str_dt.shape)
    
    #########################
    # Get indices of storms 
    # so that ids_storms[id] gives the ids in the arrays
    # str_id, str_lon,.. belonging to that specific storm
    #########################
    uniq_ids = np.unique(str_id)
    ids_storms = get_indices_sparse(str_id)
    nrstorms = len(uniq_ids)

    #########################
    # Define result arrays
    #########################
    connTracks = dok_matrix((nrstorms,nrstorms))
    angleTracks = dok_matrix((nrstorms,nrstorms))
    drTracks  = dok_matrix((nrstorms,nrstorms))
    dtTracks = dok_matrix((nrstorms,nrstorms))

    #########################
    # Preprocess storm tracks
    #########################
    #Check which hemisphere belongs storms to
    hemstorms = np.full(nrstorms,"Undefined")
    firstdt = []
    lastdt = []

    for strid in range(nrstorms):    
        dt_temp = str_dt[ids_storms[uniq_ids[strid]]]
        lat_temp = str_lat[ids_storms[uniq_ids[strid]]]

        #Save the first and last dt
        firstdt.append(dt_temp[0])
        lastdt.append(dt_temp[-1])

        #Check if the storm is in the NH or SH
        if(np.nanmean(lat_temp) > 0):
            hemstorms[strid] = "NH"
        elif(np.nanmean(lat_temp) < 0):
            hemstorms[strid] = "SH"

    firstdt = np.array(firstdt)
    lastdt = np.array(lastdt)

    # START CALCULATION OF CLUSTERS
    print("---------------------------------------------")
    print("Start checking for:                          ")
    print("Distance threshold = " + str(Options["distthresh"]))
    print("Time threshold = " + str(Options["timthresh"]))
    print("Length threshold = " + str(Options["lngthresh"]))
    print("---------------------------------------------")

    #Convert timthresh to td object 
    timthresh_dt = td(hours=Options["timthresh"])

    ######################################################
    # Step 1 Find connected and clustered storms
    #######################################################
    for strm1 in range(nrstorms): 
        if(strm1%100 == 0):
            print(strm1) 
        selidxs1 = ids_storms[uniq_ids[strm1]] 

        lats1 = str_lat[selidxs1]	
        lons1 = str_lon[selidxs1]
        times1 = str_dt[selidxs1]

        #Only compare with storms which are close enought im time compared to strm1 
        diffdt1  = firstdt - np.array(lastdt[strm1])
        diffdt2  = np.array(firstdt[strm1]) - lastdt

        #To do: Check if this can be speed up
        strm2idxs = np.where((np.arange(nrstorms) > strm1) & ((diffdt1 <= timthresh_dt) & (diffdt2 <= timthresh_dt)) & (hemstorms == hemstorms[strm1]))[0]

        for strm2 in strm2idxs: 

            selidxs2 = ids_storms[uniq_ids[strm2]] 
            lats2 = str_lat[selidxs2]
            lons2 = str_lon[selidxs2] 
            times2 = str_dt[selidxs2]

            #Check if storm 1 and 2 are connected
            conn, angle, dt, dr, strConn1, strConn2  =\
                connect_cyclones(lons1,lats1,times1,lons2,lats2,times2,Options)

            #Save Results in arrays
            connTracks[strm2,strm1] = conn
            connTracks[strm1,strm2] = conn
            angleTracks[strm1,strm2] = angle
            dtTracks[strm1,strm2] = dt
            drTracks[strm1,strm2] = dr

            str_connected[selidxs1] += strConn1
            str_connected[selidxs2] += strConn2

    #Reformat sparse matrix for efficient row slicing
    connTracks = connTracks.tocsr()

    ########################
    # Step 2 Find clusters
    ########################
    clusters = []
    maxlength = 1

    for stridx in range(nrstorms):
        #print(stridx)
        clusttemp = find_cluster_type_dokm([stridx],connTracks)        

        if(len(clusttemp) > maxlength):
            maxlength = len(clusttemp)

        clusttemp = [uniq_ids[x] for x in clusttemp] #Convert indices to storm id
        clusters.append(clusttemp)

    #Delete duplicates and sort on the first number in clusters:
    unique_clusters = [list(x) for x in set(tuple(x) for x in clusters)]

    #from operator import itemgetter
    sorted_clusters =  sorted(unique_clusters)

    ############################
    # Step 3 Suborder clusters
    ############################
    sorted_subclusters_bjerknes = []
    sorted_subclusters_stagnant = []

    for cluster in sorted_clusters:
        subclusters_bjerknes = []
        subclusters_stagnant = []

        for strid in cluster:

            #Convert strid to index
            stridx = [i for i in range(len(uniq_ids)) if uniq_ids[i] == strid]
            #np.where(uniq_ids == strid)[0]

            #Length clusters
            clusttemp = find_cluster_type_dokm(stridx,connTracks,contype="Bjerknes")

            clusttemp = [uniq_ids[x] for x in clusttemp] #Convert indices to storm id
            subclusters_bjerknes.append(clusttemp)

            #Stationary clusters
            clusttemp = find_cluster_type_dokm(stridx,connTracks,contype="Stagnant")

            clusttemp = [uniq_ids[x] for x in clusttemp] #Convert indices to storm id
            subclusters_stagnant.append(clusttemp)

        #Delete duplicates and sort on the first number in (sub)clusters:
        unique_subclusters = [list(x) for x in set(tuple(x) for x in subclusters_bjerknes)]
        sorted_subclusters_bjerknes.append(sorted(unique_subclusters))

        #Delete duplicates and sort on the first number in (sub)clusters:
        unique_subclusters = [list(x) for x in set(tuple(x) for x in subclusters_stagnant)]
        sorted_subclusters_stagnant.append(sorted(unique_subclusters))

    sorted_clusters_bjerknes = sorted(unnest(sorted_subclusters_bjerknes))
    sorted_clusters_stagnant = sorted(unnest(sorted_subclusters_stagnant))

    # return results
    return sorted_clusters, sorted_clusters_bjerknes, sorted_clusters_stagnant, str_connected

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
################################################## Required functions #######################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

##### Temporary until upgraded 
#import sys
#sys.path.append('/g/data/gb02/mb0427/wxlib_packages/lib/python3.10/site-packages/') 
##from dynlib import derivatives, gridlib#, utils, proj
#import derivatives_python as derivatives

# Configuration options
#lmsl_thres_closed = 0.2 * 1.0e-7/8.1 # conversion from hPa / deg-lat^2 to Pa/m^2 
#lmsl_thres_open = 0.6 * 1.0e-7/8.1 
#msl_min_prominence_closed = 150         # Local prominence in Pa compared to environment, as defined by filter_size and thus maxdist_lmsl_center
#msl_min_prominence_open = 25
#mindist_centers = 750.0e3
#maxdist_lmsl_center = 500.0e3   # maximmum distance between max-Laplace and either minval (closed system) or mingrad (open system)
#maxoro = 1000
#
#ftopeq = 0 #0.5                    # Reduce lmsl by loro (suitable for msl in Pa / and oro in m)
#
#maxdist_track = 500.0e3
#dp_per_dist = 300/100.0e3       # 100 km distance are equivalent to 3 hPa


def make_polar_grid(hemis, nxy=2000, transform=ccrs.PlateCarree()):
    if hemis == 'nh':

        # Set up the North-Polar Stereographic projection
        projection = ccrs.NorthPolarStereo()
        
        # Transform points to the projection coordinates
        x_nb, y_nb = projection.transform_point(180.0, 20.0, transform)  # Northern boundary
        x_sb, y_sb = projection.transform_point(0.0, 20.0, transform)  # Southern boundary
        x_np, y_np = projection.transform_point(0.0, 90.0, transform)  # North Pole
    
    else:
        # Set up the North-Polar Stereographic projection
        projection = ccrs.SouthPolarStereo()
        
        # Transform points to the projection coordinates
        x_nb, y_nb = projection.transform_point(0, -20, transform)  # Northern boundary
        x_sb, y_sb = projection.transform_point(180.0, -20.0, transform)  # Southern boundary
        x_np, y_np = projection.transform_point(0.0, -90.0, transform)  # North Pole

    
    # Average grid spacing
    dxy = (y_nb - y_sb)/nxy
    
    # Create grid in Cartesian coordinates
    xy_vals = np.arange(dxy/2 + y_sb, y_nb, dxy)
    y, x = np.meshgrid(xy_vals, xy_vals)

    # Inverse project the x, y coordinates back to lat/lon
    lonlat = transform.transform_points(projection, x, y)
    lon, lat = lonlat[..., 0], lonlat[..., 1]
    
    # For clon/clat (constant longitude at x_np, varying latitude at xy_vals)
    clonlat = transform.transform_points(projection, np.ones(xy_vals.shape) * x_np, xy_vals)
    clon, clat = clonlat[..., 0], clonlat[..., 1]
    
    # Derive actual dx, dy from the spacing of latitudes along the 0/180-meridian
    dlon = np.abs(clon[1:] - clon[:-1])
    over_np = np.where(dlon > 170.0)[0]

    dxy_real = np.abs(clat[1:] - clat[:-1])
    dxy_real[over_np] = 180.0 - np.abs(clat[over_np]) - np.abs(clat[over_np])
    dxy_real = (dxy_real[1:] + dxy_real[:-1]) * 111111.111

    dx = np.empty(lon.shape)
    dy = np.empty(lon.shape)

    dx[:,1:-1] = dxy_real[np.newaxis,:]
    dx[:,0] = dxy_real[0]
    dx[:,-1] = dxy_real[-1]

    dy[1:-1,:] = dxy_real[:,np.newaxis]
    dy[0,:] = dxy_real[0]
    dy[-1,:] = dxy_real[-1]

    grid = gridlib.grid_by_xy(x, y, dx=dx, dy=dy, lats=lat, lons=lon)

    return projection, grid 


def sort_cycpos(x1, y1, x2, y2, p1=None, p2=None, dist_thres=250.0e3, dp_per_dist=0.003):
    # Check which positions in match between the two sets of coordinates, 
    # return boolean mask of matches for both first and second input

    dists = ((x1[np.newaxis,:]-x2[:,np.newaxis])**2 + (y1[np.newaxis,:]-y2[:,np.newaxis])**2)
    
    # Optional: take into account pressure differences in the distance metric
    if not type(p1) == type(None):
        dists += (p1[np.newaxis,:]-p2[:,np.newaxis])**2 / dp_per_dist**2
        
    mindists1 = dists.min(axis=0)
    mindists2 = dists.min(axis=1)
    
    match1 = []
    match2 = []
    associations = {}
    
    for n in range(len(x1)):
        for m in range(len(x2)):
            if m in match2:
                continue

            if dists[m,n] <= dist_thres**2 and dists[m,n] == mindists1[n] and dists[m,n] == mindists2[m]:
                match1.append(n)
                match2.append(m)
                associations[n] = m
                break

    return match1, match2, associations


def _filter_locs_by_mask(locs, mask):
    outmask = mask[locs[0],locs[1]]

    return tuple(locs[i][~outmask] for i in range(len(locs)))


def _filter_cyclones_by_mindist(cycs, mindist_thres):
    lens = len(cycs)
    mask = np.ones((lens,), dtype='bool')
    for n in range(lens):
        cyc1 = cycs.iloc[n]
        for m in range(n+1,lens):
            cyc2 = cycs.iloc[m]
            dist = np.sqrt((cyc1.x - cyc2.x)**2 + (cyc1.y - cyc2.y)**2)
            if dist < mindist_thres:
                if cyc1.msl_maxlap >= cyc2.msl_maxlap:
                    mask[m] = False
                else:
                    mask[n] = False
    
    return cycs[mask]


def locate_cyclones(date, msl, grid, hemis, 
                        lmsl_thres_closed, lmsl_thres_open, maxdist_lmsl_center,
                        msl_min_prominence_closed, msl_min_prominence_open, mindist_centers,
                        ftopeq=0,maxoro=1000,**kwargs):
    # Derived options
    lmsl_thres = min(lmsl_thres_closed, lmsl_thres_open)
    dx_median = np.median(grid.dx)
    filter_size = int(maxdist_lmsl_center / dx_median * np.sqrt(2))   # Reflecting max distance between Max-Lap and these features

    
    # Look for minima/maxima
    mslx, msly = derivatives.grad(msl[np.newaxis,:,:], grid.dx, grid.dy)
    gmsl = np.sqrt(mslx[0,:,:]**2 + msly[0,:,:]**2)
    del mslx, msly
    lmsl = derivatives.lap2(msl[np.newaxis,:,:], grid.dx, grid.dy)[0,:,:]
    msl_prominence = ndimg.uniform_filter(msl, size=int(filter_size*np.sqrt(2))) - msl

    # Optional: reduce lmsl by loro by Laplacian of orography
    loro = derivatives.lap2(grid.oro[np.newaxis,:,:], grid.dx, grid.dy)[0,:,:]
    lmsl -= ftopeq * np.abs(loro)

    maxlap = np.where(np.logical_and(
        ndimg.maximum_filter(lmsl, size=filter_size) == lmsl,
        lmsl >= lmsl_thres
    ))
    minval = np.where(ndimg.minimum_filter(msl, size=filter_size) == msl)
    mingrad = np.where(ndimg.minimum_filter(gmsl, size=filter_size) == gmsl)


    # Remove minima/maxima over high topography
    maxlap = _filter_locs_by_mask(maxlap, grid.oro > maxoro)
    minval = _filter_locs_by_mask(minval, grid.oro > maxoro)
    mingrad = _filter_locs_by_mask(mingrad, grid.oro > maxoro)


    # Coordinates of the Max-Laplace, min-gradient and min-values in map coordinates and lat/lon coordinates
    x_maxlap, y_maxlap = grid.x[maxlap[0],maxlap[1]], grid.y[maxlap[0],maxlap[1]]
    x_minval, y_minval = grid.x[minval[0],minval[1]], grid.y[minval[0],minval[1]]
    x_mingrad, y_mingrad = grid.x[mingrad[0],mingrad[1]], grid.y[mingrad[0],mingrad[1]]

    lon_maxlap, lat_maxlap = grid.lon[maxlap[0],maxlap[1]], grid.lat[maxlap[0],maxlap[1]]
    lon_minval, lat_minval = grid.lon[minval[0],minval[1]], grid.lat[minval[0],minval[1]]
    lon_mingrad, lat_mingrad = grid.lon[mingrad[0],mingrad[1]], grid.lat[mingrad[0],mingrad[1]]


    # Match Max-Laplace with min-values -> closed systems
    match_minval, match_maxlap, ___ = sort_cycpos(x_minval, y_minval, x_maxlap, y_maxlap, dist_thres=maxdist_lmsl_center)
    
    # Properties of closed systems
    shape = (len(match_minval), )
    msl_center = msl[minval[0],minval[1]][match_minval]
    cyclones_closed = pd.DataFrame({
        'track_id': -np.ones(shape, dtype='i4'),                                    # Cyclone Track ID (to be assigned later)
        'hemis': [hemis,]*shape[0],                                                 # Hemisphere, NH or SH
        'date': [date,]*shape[0],                                                   # Date
        'x': x_minval[match_minval], 'y': y_minval[match_minval],                   # Location of the MSL minimum in map coordinates
        'u': np.empty(shape)*np.nan, 'v': np.empty(shape)*np.nan,                   # Movement vector in map coordinates
        'segment_len': np.empty(shape)*np.nan,                                      # Distance to previous point of track
        'movement_speed': np.empty(shape)*np.nan,                                   # Movement speed since previous point of track
        'movement_direction': np.empty(shape)*np.nan,                               # Direction from previous point of track
        'x_pred': [np.empty((2,1)) * np.nan for i in range(shape[0])],              # Kalman state vector before the update
        'y_pred': [np.empty((2,1)) * np.nan for i in range(shape[0])],              # Kalman state vector before the update
        'x_filter': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Kalman state vector
        'y_filter': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Kalman state vector
        'x_smooth': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Smoothed state vector
        'y_smooth': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Smoothed state vector
        'Px_pred': [np.empty((2,2)) * np.nan for i in range(shape[0])],             # Kalman uncertainty covaiance matrix before the update
        'Py_pred': [np.empty((2,2)) * np.nan for i in range(shape[0])],             # Kalman uncertainty covaiance matrix before the update
        'Px_filter': [np.empty((2,2)) * np.nan for i in range(shape[0])],           # Kalman uncertainty covaiance matrix
        'Py_filter': [np.empty((2,2)) * np.nan for i in range(shape[0])],           # Kalman uncertainty covaiance matrix
        'lon': np.empty(shape)*np.nan, 'lat': np.empty(shape)*np.nan,               # Location of the MSL minimum
        'lon_raw': lon_minval[match_minval], 'lat_raw': lat_minval[match_minval],   # Location of the MSL minimum
        'lon_maxlap': lon_maxlap[match_maxlap], 'lat_maxlap': lat_maxlap[match_maxlap], # Location of the Lap(MSL) maximum
        'closed': np.ones(shape, dtype='bool'),                                     # Closed cyclones are closed
        'msl': msl_center,                                                          # MSL at the location of the minimmum
        'msl_prominence': msl_prominence[minval[0],minval[1]][match_minval],        # Local prominence of the minimum MSL
        'msl_mingrad': np.zeros(shape),                                             # Min-Grad zero by definition for closed systems
        'msl_maxlap': lmsl[maxlap[0],maxlap[1]][match_maxlap],                      # Lap(MSL) at the location of the maximmum
    })


    # Match remaining Max-Laplace locations with minimum-gradients -> open systems
    remain_maxlap = np.array([i not in match_maxlap for i in range(len(x_maxlap))])
    x_remain, y_remain = x_maxlap[remain_maxlap], y_maxlap[remain_maxlap]
    match_mingrad, match_maxlap2, ___ = sort_cycpos(x_mingrad, y_mingrad, x_remain, y_remain, dist_thres=maxdist_lmsl_center)

    shape = (len(match_mingrad), )
    msl_center = msl[mingrad[0],mingrad[1]][match_mingrad]
    cyclones_open = pd.DataFrame({
        'track_id': -np.ones(shape, dtype='i4'),
        'hemis': [hemis,]*shape[0],                                                 # Hemisphere, NH or SH
        'date': [date,]*shape[0],                                                   # Date
        'x': x_mingrad[match_mingrad], 'y': y_mingrad[match_mingrad],               # Location of the min gradient in map coordinates
        'u': np.empty(shape)*np.nan, 'v': np.empty(shape)*np.nan,                   # Movement vector in map coordinates
        'segment_len': np.empty(shape)*np.nan,                                      # Distance to previous point of track
        'movement_speed': np.empty(shape)*np.nan,                                   # Movement speed since previous point of track
        'movement_direction': np.empty(shape)*np.nan,                               # Direction from previous point of track
        'x_pred': [np.empty((2,1)) * np.nan for i in range(shape[0])],              # Kalman state vector before the update
        'y_pred': [np.empty((2,1)) * np.nan for i in range(shape[0])],              # Kalman state vector before the update
        'x_filter': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Kalman state vector
        'y_filter': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Kalman state vector
        'x_smooth': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Smoothed state vector
        'y_smooth': [np.empty((2,1)) * np.nan for i in range(shape[0])],            # Smoothed state vector
        'Px_pred': [np.empty((2,2)) * np.nan for i in range(shape[0])],             # Kalman uncertainty covaiance matrix before the update
        'Py_pred': [np.empty((2,2)) * np.nan for i in range(shape[0])],             # Kalman uncertainty covaiance matrix before the update
        'Px_filter': [np.empty((2,2)) * np.nan for i in range(shape[0])],           # Kalman uncertainty covaiance matrix
        'Py_filter': [np.empty((2,2)) * np.nan for i in range(shape[0])],           # Kalman uncertainty covaiance matrix
        'lon': np.empty(shape)*np.nan, 'lat': np.empty(shape)*np.nan,               # Location of the MSL minimum
        'lon_raw': lon_mingrad[match_mingrad], 'lat_raw': lat_mingrad[match_mingrad], # Location of the MSL minimum
        'lon_maxlap': lon_maxlap[remain_maxlap][match_maxlap2], 'lat_maxlap': lat_maxlap[remain_maxlap][match_maxlap2],
        'closed': np.zeros(shape, dtype='bool'),
        'msl': msl_center,                                                          # MSL the location of the min gradient
        'msl_prominence': msl_prominence[mingrad[0],mingrad[1]][match_mingrad],     # Local prominence of the minimum MSL
        'msl_mingrad': gmsl[mingrad[0],mingrad[1]][match_mingrad],                  # Min gradient at the location of the min-gradient
        'msl_maxlap': lmsl[maxlap[0],maxlap[1]][remain_maxlap][match_maxlap2],
    })
    
    # Filter closed and open systems by additional Lap(MSL) and other thresholds
    below_minprom = cyclones_closed.msl_prominence <= msl_min_prominence_closed
    closed2open = cyclones_closed[below_minprom]
    closedmask = np.logical_and(cyclones_closed.msl_maxlap >= lmsl_thres_closed, ~below_minprom)
    cyclones_closed = cyclones_closed[closedmask]

    openmask = np.logical_and(
            cyclones_open.msl_maxlap >= lmsl_thres_open, 
            cyclones_open.msl_prominence >= msl_min_prominence_open
    )
    cyclones_open = cyclones_open[openmask]

    openmask2 = np.logical_and(
            closed2open.msl_maxlap >= lmsl_thres_open, 
            closed2open.msl_prominence >= msl_min_prominence_open
    )
    closed2open = closed2open[openmask2]
    
    # Remove minor systems by minimum distance threshold
    cyclones = pd.concat([cyclones_closed, closed2open, cyclones_open], ignore_index=True)
    cyclones = _filter_cyclones_by_mindist(cyclones, mindist_centers)
    cyclones = cyclones.reset_index(drop=True)
    
    return cyclones


# Physcial model used in the Kalman filter and RTS smoother: constant velocity
tstep = 10800                      # Time step of 3h
F = np.array([[1, tstep], [0, 1]]) # Operator
R = (90.0e3)**2                 # Measurement uncertainty constant at 90 km
Q = np.array([[tstep**2, tstep], [tstep, 1]]) * 25 # Process noise error of 5 m/s, i.e. how bad is the assumption of constant-velocity?
H = np.array([[1, 0],])         # Only position is observed

def kalman_step(x, Px, zx):
    # 1a. Extrapolate state
    x_ = F @ x
    # 1b. Extrapolate uncertainty
    Px_ = F @ Px @ F.T + Q

    # 2a. Kalman gain
    K = Px_ @ H.T / (H @ Px_ @ H.T + R)
    # 2b. Update state
    x = x_ + K * (zx - H @ x_)
    # 2c. Update uncertainty
    ImKH = np.eye(Px_.shape[0]) - K @ H
    Px = ImKH @ Px_ @ ImKH.T + R * K @ K.T

    return x, Px, x_, Px_, K


# Rauch-Tung-Striebel smoother
def rts_smoother_step(x, Px, x_, Px_, xs, Pxs):
    C = Px @ F.T @ np.linalg.inv(Px_)
    
    xs = x + C @ (xs - x_)
    Pxs = Px + C @ (Pxs - Px_) @ C.T

    return xs, Pxs, C


def smooth_cyclone_track(maps, all_cyc, cyc, tstep, 
                         transform=ccrs.PlateCarree()):
    # Initialisation
    hemis = cyc.iloc[-1].hemis
    xs = cyc.iloc[-1].x_filter
    Pxs = cyc.iloc[-1].Px_filter
    ys = cyc.iloc[-1].y_filter
    Pys = cyc.iloc[-1].Py_filter

    #lon, lat = maps[hemis](xs[0,0], ys[0,0], inverse=True)
    #print(xs)
    #lon, lat = transform.transform_point(xs[0,0], ys[0,0], maps[hemis])
    lonlat = transform.transform_points(maps[hemis], xs[0,0], ys[0,0])
    lon, lat = lonlat[..., 0], lonlat[..., 1]
    
    _save_to_cyc(all_cyc, cyc.index[-1], **dict(x_smooth=xs, y_smooth=ys, lon=lon, lat=lat))
    prevlon, prevlat = lon, lat

    for i in range(len(cyc)-2,-1,-1):
        xs, Pxs, C = rts_smoother_step(
                cyc.iloc[i].x_filter, cyc.iloc[i].Px_filter, 
                cyc.iloc[i+1].x_pred, cyc.iloc[i+1].Px_pred,
                xs, Pxs)
        ys, Pys, C = rts_smoother_step(
                cyc.iloc[i].y_filter, cyc.iloc[i].Py_filter, 
                cyc.iloc[i+1].y_pred, cyc.iloc[i+1].Py_pred,
                ys, Pys)
        
        #lon, lat = maps[hemis](xs[0,0], ys[0,0], inverse=True)
        lon, lat = transform.transform_point(xs[0,0], ys[0,0], maps[hemis])
        _save_to_cyc(all_cyc, cyc.index[i], **dict(x_smooth=xs, y_smooth=ys, lon=lon, lat=lat))

        segment_len = utils.dist_sphere(lon, lat, prevlon, prevlat)
        speed = segment_len / tstep
        bearing = utils.direction_on_sphere(lon, lat, prevlon, prevlat)

        _save_to_cyc(all_cyc, cyc.index[i+1], **dict(segment_len=segment_len, movement_speed=speed, movement_direction=bearing))
        prevlon, prevlat = lon, lat
    
    return

def _save_to_cyc(cyc, cycid, **kwargs):
    for key, value in kwargs.items():
        if type(value) == np.ndarray:
            if np.isscalar(cyc.loc[cycid, key]):
                cyc.loc[cycid, key] = value  # Direct assignment for scalar values
            else:
                cyc.loc[cycid,key][:,:] = value
        else:
            cyc.loc[cycid,key] = value
            
    return

def track_cyclones(ntracks, cyclones_prev, cyclones_cur, maxdist_track=750.0e3, **kwargs):
    # Match previous with current cyclone positions
    x_prev, y_prev, msl_prev = cyclones_prev.x.to_numpy(), cyclones_prev.y.to_numpy(), cyclones_prev.msl.to_numpy()
    u_prev, v_prev = cyclones_prev.u.to_numpy(), cyclones_prev.v.to_numpy()
    x_cur, y_cur, msl_cur = cyclones_cur.x.to_numpy(), cyclones_cur.y.to_numpy(), cyclones_cur.msl.to_numpy()
    
    # Predict cyclone positions based on prior movement
    notnan = ~np.isnan(u_prev)
    x_prev[notnan] += u_prev[notnan] * tstep
    y_prev[notnan] += v_prev[notnan] * tstep

    ___, ___, associations = sort_cycpos(x_prev, y_prev, x_cur, y_cur, msl_prev, msl_cur, 
            dist_thres=maxdist_track, dp_per_dist = 300/100.0e3)
    
    for id_prev, id_cur in associations.items():
        if cyclones_prev.loc[id_prev,'track_id'] >= 0:
            # Third-to-final Kalman step
            x, Px, x_, Px_, Kx = kalman_step(
                            cyclones_prev.loc[id_prev,'x_filter'], 
                            cyclones_prev.loc[id_prev,'Px_filter'], 
                            cyclones_cur.loc[id_cur,'x'])
            y, Py, y_, Py_, Ky = kalman_step(
                            cyclones_prev.loc[id_prev,'y_filter'], 
                            cyclones_prev.loc[id_prev,'Py_filter'], 
                            cyclones_cur.loc[id_cur,'y'])
            
            _save_to_cyc(cyclones_cur, id_cur, **dict(
                    track_id=cyclones_prev.loc[id_prev,'track_id'], u=x[1,0], v=y[1,0], 
                    x_filter=x, y_filter=y, Px_filter=Px, Py_filter=Py,
                    x_pred=x_, y_pred=y_, Px_pred=Px_, Py_pred=Py_)
            )


        else:
            uraw = (cyclones_cur.loc[id_cur,'x'] - cyclones_prev.loc[id_prev,'x'])/tstep
            vraw = (cyclones_cur.loc[id_cur,'y'] - cyclones_prev.loc[id_prev,'y'])/tstep

            # Initialise Kalman filter: rest at origin, but immense uncertainty in x, y, u, v
            Px = np.array([[1.0e6, 0], [0, 100]])
            Py = np.array([[9.9e99, 0], [0, 9.9e99]])
            
            # First Kalman step
            x = np.array([[cyclones_prev.loc[id_prev,'x']], [uraw]])
            y = np.array([[cyclones_prev.loc[id_prev,'y']], [vraw]])
            x, Px, x_, Px_, Kx = kalman_step(x, Px, cyclones_prev.loc[id_prev,'x'])
            y, Py, y_, Py_, Ky = kalman_step(y, Py, cyclones_prev.loc[id_prev,'y'])

            _save_to_cyc(cyclones_prev, id_prev, **dict(
                    track_id=ntracks, u=uraw, v=vraw,
                    x_filter=x, y_filter=y, Px_filter=Px, Py_filter=Py,
                    x_pred=x_, y_pred=y_, Px_pred=Px_, Py_pred=Py_)
            )

            # Second Kalman step
            x, Px, x_, Px_, Kx = kalman_step(x, Px, cyclones_cur.loc[id_cur,'x'])
            y, Py, y_, Py_, Ky = kalman_step(y, Py, cyclones_cur.loc[id_cur,'y'])

            _save_to_cyc(cyclones_cur, id_cur, **dict(
                    track_id=ntracks, u=uraw, v=vraw, 
                    x_filter=x, y_filter=y, Px_filter=Px, Py_filter=Py,
                    x_pred=x_, y_pred=y_, Px_pred=Px_, Py_pred=Py_)
            )

            ntracks += 1

    return ntracks


def aggregate_trackinfo(cyc_locs):
    # Derive some track diagnostics from cyclone locations
    genesis = cyc_locs.iloc[0].date
    lysis = cyc_locs.iloc[-1].date
    lifetime = (lysis - genesis).total_seconds()/3600.0

    x = np.array([val[0,0] for val in cyc_locs.x_smooth])
    y = np.array([val[0,0] for val in cyc_locs.y_smooth])
    dist = np.sqrt((x[-1]-x[0])**2 + (y[-1]-y[0])**2)
    segment_lens = [np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2) for i in range(1,len(x))]

    maxlap = cyc_locs[cyc_locs.msl_maxlap == cyc_locs.msl_maxlap.max()].iloc[0]
    minmsl = cyc_locs[cyc_locs.msl == cyc_locs.msl.min()].iloc[0]
    maxprom = cyc_locs[cyc_locs.msl_prominence == cyc_locs.msl_prominence.max()].iloc[0]
        
    # Collect information in (single-row) dataframe (to be concatenated later)
    trackinfo = pd.DataFrame({
        'track_id': cyc_locs.iloc[0].track_id,              # Cyclone Track ID
        'lifetime': lifetime,                               # Cyclone life time [hours]
        'distance_firstlast': dist,                         # Distance covered (between last and first point) [km]
        'distance_traveled': sum(segment_lens),             # Length of track over all points, i.e. total distance travelled
        'closed': np.any(cyc_locs.closed),                  # Did cyclone mature to closed status?

        'date_genesis': genesis,                            # Genesis date and location
        'lon_genesis': cyc_locs.iloc[0].lon, 'lat_genesis': cyc_locs.iloc[0].lat,
        'date_lysis': lysis,                                # Lysis date and location
        'lon_lysis': cyc_locs.iloc[-1].lon, 'lat_lysis': cyc_locs.iloc[-1].lat,

        'msl_maxlap': maxlap.msl_maxlap,                    # Value of maximum Laplacian [Pa/km^2]
        'date_maxlap': maxlap.date,                         # Date and location of maximum Laplacian
        'lon_maxlap': maxlap.lon, 'lat_maxlap': maxlap.lat, 
        'msl_min': minmsl.msl,                              # Value of minimum MSL [Pa]
        'date_mslmin': minmsl.date,                         # Date and location of minimum MSL
        'lon_mslmin': minmsl.lon, 'lat_mslmin': minmsl.lat, 
        'msl_maxprom': maxprom.msl_prominence,              # Value of maximum local prominence [Pa]
        'date_maxprom': maxprom.date,                       # Date and location of maximum local prominence 
        'lon_maxprom': maxprom.lon, 'lat_maxprom': maxprom.lat,
    }, index=[0,])                                          # Index will be reset on concatenation

    return trackinfo


def filter_tracks(tracks, lastdate, **kwargs):
    # Filter tracks by: 
    # - minimum lifetime of 24 hours, 
    # - minimum distance between first and last point of 500 km, and 
    # - must be closed at some point
    filter_fun = lambda trk: (
        trk.lifetime >= 24 and
        trk.distance_firstlast >= 500.0e3 and
        trk.closed == True
    )

    unfinished = np.array([trk.date_lysis == lastdate for trkid, trk in tracks.iterrows()])

    mask = np.array([filter_fun(trk) for trkid, trk in tracks.iterrows()])
    mask[unfinished] = False

    return tracks.iloc[mask], tracks.iloc[unfinished]