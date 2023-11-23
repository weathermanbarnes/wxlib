#!/usr/bin/env python
# -*- encoding: utf-8

import numpy as np
import pandas as pd
import scipy.interpolate as interp
import scipy.ndimage as ndimg

from . import derivatives, utils, gridlib, proj


# Configuration options
lmsl_thres_closed = 0.2 * 1.0e-7/8.1 # conversion from hPa / deg-lat^2 to Pa/m^2 
lmsl_thres_open = 0.6 * 1.0e-7/8.1 
msl_min_prominence_closed = 150         # Local prominence in Pa compared to environment, as defined by filter_size and thus maxdist_lmsl_center
msl_min_prominence_open = 25
mindist_centers = 750.0e3
maxdist_lmsl_center = 500.0e3   # maximmum distance between max-Laplace and either minval (closed system) or mingrad (open system)
maxoro = 1000

ftopeq = 0 #0.5                    # Reduce lmsl by loro (suitable for msl in Pa / and oro in m)

maxdist_track = 500.0e3
dp_per_dist = 300/100.0e3       # 100 km distance are equivalent to 3 hPa


def make_polar_grid(hemis, nxy=2000):
    if hemis == 'nh':
        m = proj.Basemap(projection='npstere', boundinglat=20, lon_0=0)

        # "northern" and "southern" boundaries as well as North Pole in Cartiesian projection coordinates
        x_nb, y_nb = m(180.0, 20.0)
        x_sb, y_sb = m(0.0, 20.0)
        x_np, y_np = m(0.0, 90.0)
    else:
        m = proj.Basemap(projection='spstere', boundinglat=-20, lon_0=0)
    
        # "northern" and "southern" boundaries as well as North Pole in Cartiesian projection coordinates
        x_nb, y_nb = m(0.0, -20.0)
        x_sb, y_sb = m(180.0, -20.0)
        x_np, y_np = m(0.0, -90.0)

    
    # Average grid spacing
    dxy = (y_nb - y_sb)/nxy
    
    # Create grid in Cartesian coordinates
    xy_vals = np.arange(dxy/2 + y_sb, y_nb, dxy)
    y, x = np.meshgrid(xy_vals, xy_vals)
    
    # Inverse projection to lat/lon coordinates
    lon, lat = m(x, y, inverse=True)
    clon, clat = m(np.ones(xy_vals.shape)*x_np, xy_vals, inverse=True)
    
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

    return m, grid 


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


def locate_cyclones(date, msl, grid, hemis):
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


def smooth_cyclone_track(maps, all_cyc, cyc):
    # Initialisation
    hemis = cyc.iloc[-1].hemis
    xs = cyc.iloc[-1].x_filter
    Pxs = cyc.iloc[-1].Px_filter
    ys = cyc.iloc[-1].y_filter
    Pys = cyc.iloc[-1].Py_filter

    lon, lat = maps[hemis](xs[0,0], ys[0,0], inverse=True)
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
        
        lon, lat = maps[hemis](xs[0,0], ys[0,0], inverse=True)
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
            cyc.loc[cycid,key][:,:] = value
        else:
            cyc.loc[cycid,key] = value
    
    return


def track_cyclones(ntracks, cyclones_prev, cyclones_cur):
    # Match previous with current cyclone positions
    x_prev, y_prev, msl_prev = cyclones_prev.x.to_numpy(), cyclones_prev.y.to_numpy(), cyclones_prev.msl.to_numpy()
    u_prev, v_prev = cyclones_prev.u.to_numpy(), cyclones_prev.v.to_numpy()
    x_cur, y_cur, msl_cur = cyclones_cur.x.to_numpy(), cyclones_cur.y.to_numpy(), cyclones_cur.msl.to_numpy()
    
    # Predict cyclone positions based on prior movement
    notnan = ~np.isnan(u_prev)
    x_prev[notnan] += u_prev[notnan] * tstep
    y_prev[notnan] += v_prev[notnan] * tstep

    ___, ___, associations = sort_cycpos(x_prev, y_prev, x_cur, y_cur, msl_prev, msl_cur, 
            dist_thres=maxdist_track, dp_per_dist=dp_per_dist)
    
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


def filter_tracks(tracks, lastdate):
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


# C'est le fin
