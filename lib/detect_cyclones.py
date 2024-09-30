##### Import packages

import sys
import pickle
import numpy as np
import scipy.interpolate as interp
from datetime import datetime as dt, timedelta as td

##### Import created packages
import utils

#####

def cyclone_by_lapmsl(msl, grid, prev_cyc=None, prev_tracks=None, quiet=False, lmsl_thres_closed=2.0e-09, lmsl_thres_open=7.0e-09, maxdist_lmsl_center=500.0e3, msl_min_prominence_closed=150, msl_min_prominence_open=25, mindist_centers=750.0e3,maxdist_track=750.0e3,tstep=10800,min_lifetime=24,min_dist_travelled=500.0e3):
    
    ''' Detect and track cyclone centers following the Melbourne algorithm
    DEVELOPED BY: Clemens Spensberger (Bergen), Adapted for wxlib by Michael Barnes (Monash)

   The original detection algorithm is defined in Murray and Simmonds (1991a,b; MS91). This is a simplified
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
   
    Returns
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
    
    from cyclone_by_lapmsl_helpers import pd, make_polar_grid, locate_cyclones, track_cyclones, \
                smooth_cyclone_track, aggregate_trackinfo, filter_tracks
    
    # Prepare grid objects for the polar stereographic grids on which cyclones are being detected and tracked
    timer = dt.now()
    ilons, ilats = utils.concat1lonlat(grid.x, grid.y)
    ilats = ilats[::-1,0]
    ilons = ilons[0,:]
    ifunc = interp.RectBivariateSpline(ilats, ilons, utils.concat1(grid.oro[::-1,:]))
    
    maps = {}
    grids = {}
    for hemis in ['nh', 'sh']:
        maps[hemis], grids[hemis] = make_polar_grid(hemis)
        grids[hemis].oro = ifunc(grids[hemis].lat, grids[hemis].lon, grid=False)
    if not quiet:
        print(f'Preparation: {(dt.now()-timer).total_seconds():6.1f} seconds.')
    
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
            print(f'Detection & tracking: {tidx:03d}/{msl.shape[0]:03d}', end=chr(13))
        ifunc = interp.RectBivariateSpline(ilats, ilons, utils.concat1(msl[tidx,::-1,:]))

        for hemis in ['nh', 'sh']:
            hemis_msl = ifunc(grids[hemis].lat, grids[hemis].lon, grid=False)

            new_cyclones = locate_cyclones(grid.t_parsed[tidx], hemis_msl, grids[hemis], hemis,\
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
        print(f'Detection & tracking: {(dt.now()-timer).total_seconds():6.1f} seconds.')

    timer = dt.now()
    # Now that the basic detection and tracking is complete: run the Kalman smoother for a best estimate
    # of cyclone position and velocity
    tracks = []
    trkids = set(cyclones.track_id)
    if not quiet:
        print('Smoothing & trackinfo.', end=chr(13))
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
        print(f'Smoothing & trackinfo: {(dt.now()-timer).total_seconds():6.1f} seconds.')

    timer = dt.now()
    # Filter tracks for minimum length, minimum movement, etc.
    if not quiet:
        print('Filtering tracks.              ', end=chr(13))
    tracks, unfinished_tracks = filter_tracks(tracks, lastdate=grid.t_parsed[-1],min_lifetime=min_lifetime, min_dist_travelled=min_dist_travelled)

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