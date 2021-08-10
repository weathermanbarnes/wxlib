#!/usr/bin/env python
# -*- encoding: utf-8

''' Utilities that don't fit elsewhere '''

from __future__ import absolute_import, unicode_literals, print_function

import sys

from . import dynfor
from . import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.utils, 'utils', sys.modules[__name__])


import math
import numpy as np
import scipy as sp
import scipy.interpolate as intp
from scipy.special import erfinv

from . import tagg

from datetime import datetime as dt, timedelta as td
import calendar



def scale(var, cut=slice(None)):
    """ Automatic scaling according to netcdf attributes ``scale_factor`` and ``add_offset``

    If present missing or fill values, taken from the netCDF attributes ``missing_value``
    and ``_FillValue`` are converted to NaN.

    Parameters
    ----------

    var : nc.NetCDFVariable
        The variable to be scaled
    cut : slice
        *Optional*, default ``slice(None)``. Limit the request to a given time slice.
        The cut is always applied to the first dimension in the dat, so make sure to
        use the keyword only with time-dependent data.
        With the default data layout, only relevant data needs to be read when only 
        a time slice of the entire data is requested. Hence, using cut to limit your 
        data request can make reading the data largely more efficient.
    
    Returns
    -------
    np.ndarray
        The scaled data contained in the nc.NetCDFVariable object `var`.
    """
    
    # Apply cut first, for speed
    #if len(var.shape) >= 2:
    #    dat = var[cut,::]
    #else:
    #    dat = var[::]
    dat = var[cut]
    dat = dat.astype('f8')

    # Mask missing and fill values
    if hasattr(var, 'missing_value'):
        dat[dat == var.missing_value] = np.nan
    if hasattr(var, '_FillValue'):
        dat[dat == var._FillValue] = np.nan

    # Apply scaling, numpy is faster than Fortran function!
    if hasattr(var, 'scale_factor') or hasattr(var, 'add_offset'):
        dat = dat*getattr(var, 'scale_factor', 1.0) + getattr(var, 'add_offset', 0.0)
    
    return dat


def unscale(var):
    """ Inverse of the :func:`scale` function. 

    Compresses floating point data to 16-bit integers plus an 64-bit offset and scaling factor.

    NaN values are converted to a missing/fill value of -32767.

    Parameters
    ----------

    var : np.ndarray 
        Data to be compressed.
    
    Returns
    -------
    np.ndarray with dtype int16
        Compressed data.
    float
        Scale factor.
    float
        Offset.
    int 
        Missing value.
    """
    
    missing = -32768
    
    # Replace infinite by NaN, otherwise the scaling will reduce all other values to zero
    var = var.astype('f8')
    var[~np.isfinite(var)] = np.nan 

    maxv  = np.nanmax(var)
    minv  = np.nanmin(var)

    # divide in 2^16-2 intervals, values from -32767 -> 32767 ; reserve -32768 as missing value
    scale = (maxv-minv)/65534.0

    # Avoid division by zero in case of a constant input field
    if scale == 0.0:
        scale = 1.0

    off   = +32767*scale + minv

    res = np.round((var[::] - off)/scale)
    # Avoid integer overflow due to roundoff error
    res[res > 32737] = 32737
    res[res <= missing] = missing + 1

    res[np.isnan(res)] = missing

    return res.astype('i2'), scale, off, missing


def concat1(data):
    ''' Concatenate one latitude band in x-direction to a data array

    To be able to plot circumpolar plots without a data gap, the sphericity of the data
    must explicitly be demonstrated by contatenating the data from lon=180E to appear 
    also as data for lon=180W.

    Parameters
    ----------
    data : np.ndarray with 1-4 dimensions
        Data to be extended
    
    Returns
    -------
    np.ndarray
        Extended data
    '''
    
    if len(data.shape) == 1:
        data = np.concatenate((data, np.reshape(data[0], (1,)) ), axis=0)
    elif len(data.shape) == 2:
        data = np.concatenate((data, np.reshape(data[:,0], (data.shape[0], 1)) ), axis=1)
    elif len(data.shape) == 3:
        data = np.concatenate((data, np.reshape(data[:,:,0], (data.shape[0], data.shape[1], 1)) ), axis=2)
    elif len(data.shape) == 4:
        data = np.concatenate((data, np.reshape(data[:,:,:,0], (data.shape[0], data.shape[1], data.shape[2], 1)) ), axis=3)
    else:
        raise NotImplementedError('Concatenation not implemented for %d dimensions' % len(data.shape))
    
    return data


def concat1lonlat(x, y):
    ''' Concatenate one latitude band in x-direction to coordinate arrays

    To be able to plot circumpolar plots without a data gap, the sphericity of the data
    must explicitly be demonstrated by contatenating the data from lon=180E to appear 
    also as data for lon=180W.

    Parameters
    ----------
    x : np.ndarray with dimensions (y,x)
        Longitudes for each grid point
    y : np.ndarray with dimensions (y,x)
        Latitudes for each grid point
    
    Returns
    -------
    2-tuple of np.ndarray
        Extended coordinate arrays
    '''
    
    # Some map projections need the lon, lat array to be C-aligned.
    lon = np.ascontiguousarray(concat1(x))
    lat = np.ascontiguousarray(concat1(y))

    lon[:,-1] += 360.0

    return lon, lat


def cal_mfv(hist, bins):
    ''' Find the most frequent value from in a given histogram

    The most frequent value is defined here as the mean of the bin interval 
    boundaries.

    Parameters
    ----------
    hist : np.ndarray with 3 dimensions
         Histogram matching ``bins``.
    bins : np.ndarray or list 
         Bin intervals used for the histogram. 
    
    Returns
    -------
    np.ndarray with 2 dimensions
        Most frequently value for each grid point.
    '''
    
    s = hist[-1,:,:].sum()
    if s > 0:
        print('WARNING: hist[-1,:,:].sum() larger than zero!')

    shape = hist.shape[1:]
    mfv = np.zeros(shape)
    for j in range(shape[0]):
        for i in range(shape[1]):
            bi = hist[:-1,j,i].argmax()
            mfv[j,i] = (bins[bi+1]+bins[bi])/2.0
    
    return mfv


# Unflatten a flattened front array, using the froff list; 
#  separately for cold/warm and stationary fronts
def unflatten_fronts(fronts, froff, minlength=1):
    ''' To be made obsolete by saving cold/warm/stat fronts separately as lines in the standard-dynlib way '''

    cold = [__unflatten_fronts_t(fronts[t,0,:,:], froff[t,0,:], minlength) for t in range(fronts.shape[0])]
    warm = [__unflatten_fronts_t(fronts[t,1,:,:], froff[t,1,:], minlength) for t in range(fronts.shape[0])]
    stat = [__unflatten_fronts_t(fronts[t,2,:,:], froff[t,2,:], minlength) for t in range(fronts.shape[0])]

    return cold, warm, stat

# unflatten one time step and front type
def __unflatten_fronts_t(fronts, froff, minlength):
    ''' To be made obsolete by saving cold/warm/stat fronts separately as lines in the standard-dynlib way '''

    fronts = [fronts[froff[i]:froff[i+1],:] for i in range(froff.shape[0]-1)]
    fronts = filter(lambda lst: len(lst) >= minlength, fronts)

    return fronts

def unflatten_lines(lines, loff, static, convert_grididx=True):
    ''' Convert dynlib line format into a list of lines 

    Line coordinates (grid point indexes) are converted actual coordinates like lat/lon.
    
    '''

    xlen = static.x.shape[1]
    ylen = static.y.shape[0]
    intpx = intp.RectBivariateSpline(range(1,xlen+1), range(1,ylen+1), static.x.T)
    intpy = intp.RectBivariateSpline(range(1,xlen+1), range(1,ylen+1), static.y.T)

    lines_ = []
    for i in range(loff.shape[0]-1):
        line = lines[loff[i]:loff[i+1],:]
        if line.shape[0] < 1:
            break
        line_ = np.empty(line.shape)
        
        if convert_grididx:
            # Interpolate x and y coordinates at (fractional) grid point indexes
            line_[:,0] = intpx.ev(line[:,0], line[:,1])
            line_[:,1] = intpy.ev(line[:,0], line[:,1])
        else:
            # No conversion or interpolation
            line_[:,0:2] = line[:,0:2]
        # Copy over additional info
        line_[:,2:] = line[:,2:]

        lines_.append(line_)

    return lines_

#
# return a 3d boolean array where frontal points are True, elsewere False
def mask_fronts(fronts, froff, shape):
    ''' To be made obsolete by saving cold/warm/stat fronts separately as lines in the standard-dynlib way 
    
    See also
    --------
    ``mask_lines``, ``smear_lines``.
    '''
    masks = [np.zeros((len(fronts), shape[0], shape[1]), dtype='bool'),
         np.zeros((len(fronts), shape[0], shape[1]), dtype='bool'),
         np.zeros((len(fronts), shape[0], shape[1]), dtype='bool') ]

    for t in range(len(fronts)):
        for f in range(3):
            for n in range(froff[t,f].max()):
                # python starts counting at zero, unlike fortran
                j = round(fronts[t,f,n,1] -1)
                i = round(fronts[t,f,n,0] -1) % shape[1]
                masks[f][t,j,i] = True

    return masks


def mask_lines_with_data(lines, loff, shape, dat=None):
    ''' Mask lines in a gridded map

    Instead of returning the value ``1`` for grid points containing a line, 
    this function sets the value to either the additional info in the 
    line data or to the value of a given data array.

    Notes
    -----
    
     * Reimplement in fortran, along with hte mask_lines function there.
     * Grid size should become a parameter in the configuration

    Parameters
    ----------
    lines : np.narray with dimensions (pointindex,infotype)
        Lines to be marked on the map
    loff : np.array with dimensions (lineindex)
        List of point indexes for the first points of each line
    shape : 2-tuple of int
        Grid dimensions (ny,nx)
    dat : np.ndarray with dimensions (y,x)
        Optional: Data to be used for marking on the map
    
    Returns
    -------
    np.ndarray
        Gridded map of lines
    '''
    
    mask = np.zeros((lines.shape[0], shape[0], shape[1]))

    for t in range(lines.shape[0]):
        for n in range(loff[t].max()):
            # python starts counting at zero, unlike fortran
            j = int(round(lines[t,n,1] -1))
            i = int(round(lines[t,n,0] -1)) % shape[1]
            if type(dat) == np.ndarray:
                mask[t,j,i] = dat[t,j,i]
            else:
                mask[t,j,i] = lines[t,n,2]

    return mask


def mk_gauss(x0,stddev):
    ''' Create a Gaussian distribution function

    Parameters
    ----------
    x0 : float
        Center of the distribution
    stddev : float
        Standard deviation of the distribution
    
    Returns
    -------
    callable
        Function evaluating the Gaussian distribution function based on the given parameters
    '''

    return lambda x: np.exp(-0.5*(x-x0)**2/stddev**2)/(np.sqrt(2*np.pi)*stddev)


def smear_lines(lines, loff, static):
    ''' Smoothed and normalized line detection frequencies on a gridded map

    Grid points containing a line will be marked with the length of the line 
    through that grid point per grid point area; otherwise zero.

    Notes
    -----
    
     * Make the filter configurable

    Parameters
    ----------
    lines : np.narray with dimensions (pointindex,infotype)
        Lines to be marked on the map
    loff : np.array with dimensions (lineindex)
        List of point indexes for the first points of each line
    static : grid object
        Grid information as returned, for example, by metopen
    
    Returns
    -------
    np.ndarray
        Gridded map of lines
    '''

    filtr_len = 5
    filtr_func = mk_gauss(0, 1)
    filtr = np.array([filtr_func(x) for x in range(-filtr_len,filtr_len+1)])
    filtr /= sum(filtr)

    mask = dynfor.utils.normalize_lines(lines, loff, static.dx, static.dy)
        
    return dynfor.utils.filter_xy(mask, filtr)


def mask_insignificant(dat, mean, sig, nsig):
    ''' Mask parts of a composite that are insignificant at the chosen level

    Parameters
    ----------
    dat : np.ndarray
        Composite field
    mean : np.ndarray
        Corresponding field of mean values
    sig : np.ndarray
        Corresponding field of standard deviation values
    nsig : scalar number
        Statistical significance level, expressed normalised to the standard-derivation
    
    Returns
    -------
    np.ndarray
        Mask where the composite is significant
    '''
    dat[np.logical_and(dat < mean + nsig*sig, dat > mean - nsig*sig)] = np.nan

    return dat


igauss = lambda p: np.sqrt(2)*erfinv(2*p-1.0)
''' Inverse of the cumulative distribution function (CDF) of the Gaussian distribution 

Parameters
----------
p : scalar number
    Probability of occurrence in the range 0 - 1.

Returns
-------
float
    Threshold value in the CDF to exceed the given probability
'''


#
# 
def sect_gen_points(coords, m, dxy):
    ''' Generate (interpolation) points along a cross section
    
    Parameters
    ----------
    coords : list of 2-tuples
        List of coordinates describing the section
    m : Basemap instance
        Map projection underlying the section
    dxy : scalar number
        Maximum distance between interpolation points along the section
    
    Returns
    -------
    list
        Longitudes of the interpolation points
    list 
        Latitudes of the interpolation points
    list
        Distances along section from the section origin
    '''

    noproj = not m or type(m) == type(np)

    retlon = []
    retlat = []
    retxy  = []
    prevxy = 0
    for i, (lon1, lat1) in zip(range(len(coords)-1), coords[:-1]):
        lon2, lat2 = coords[i+1]
        
        if noproj:
            x1, x2, y1, y2 = lon1, lon2, lat1, lat2
        else: 
            (x1,x2), (y1,y2) = m((lon1,lon2), (lat1,lat2))
        d12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        N = np.ceil(d12/dxy)

        x  = [x2*alpha + x1*(1-alpha) for alpha in np.arange(0.0,(N+1)/N,1.0/N)]
        y  = [y2*alpha + y1*(1-alpha) for alpha in np.arange(0.0,(N+1)/N,1.0/N)]
        xy = [d12*alpha + prevxy for alpha in np.arange(0.0,(N+1)/N,1.0/N)]

        prevxy = xy[-1]

        if noproj:
            lon, lat = x, y
        else:
            lon, lat = m(x, y, inverse=True)
        retlon.extend(lon)
        retlat.extend(lat)
        retxy.extend(xy)

    return retlon, retlat, retxy

def aggregate(dates, dat, agg, epoch=None, bins=None):
    ''' General temporal aggregation of data

    The function assumes that the data is be equally spaced in time. This
    requirement might be relaxed in the future.

    Parameters
    ----------
    dates : list of datetime
        Dates at with the data is defined
    dat : np.ndarray
        Data to be aggregated
    agg : str
        Identifier for the aggregation period. Currently the following identifers are
        implemented: 

        * ``all``: Temporal average everything
        * ``cal_year``: Years following the calendar definition
        * ``met_season``: Seasons after their standard meteorological definition
        * ``cal_month``: Months following the calendar definition
        * ``cal_week``: Weeks following the calendar definition, starting on Monday.
        * ``cal_pendad``: Pentads according to their definition

        In addition, the following regular periodic intervals are defined:

        * ``ten_daily``: 10-day intervals starting from the first given date
        * ``weekly``: 7-day intervals starting from the first given date, 
          irrespective of the day of the week
        * ``five_daily``: 5-day intervals starting from the first given date
        * ``three_daily``: 3-day intervals starting from the first given date
        * ``two_daily``: 2-day intervals starting from the first given date
        * ``daily``: 1-day intervals starting from the first given date
    epoch : datetime
        *Optional*: Reference date used for some aggregation intervals. Defaults to the 
        first date in dates.
    bins : list of float
        *Optional*. If provided, the data will not be averaged in time. Instead histograms and 
        relatively most frequent value will be determined for each time interval.
    
    Returns
    -------
    list of datetime
        Dates at with the aggregation periods start
    list of datetime
        Dates at with the aggregation periods end
    np.ndarray
        Aggregated data averages, or most frequent values for binned data
    np.ndarray
        Number of valid data, or data histograms for binned data
    np.nparray
        List of number of timesteps included in each aggregation period
    '''

    dtd = dates[1] - dates[0]

    t_iter = getattr(tagg, agg)(dates[0], dtd, epoch=epoch)
    
    # 1. Find length of output array and generate time slices for each aggregation interval
    #    Aggegate by output of functions t_iter.start and t_iter.end; 
    #    These functions give the first and the last time step for the interval `date` belongs to
    previ = -1
    start_dates = []
    end_dates = []
    cnt_out = []
    tslc = []
    for i, date in zip(range(len(dates)), dates):
        # Previous interval unfinished, yet we are in a new one.
        # -> We apparently jumped over a range of dates and have to find a new start
        if previ >= 0 and not t_iter.start(date) == t_iter.start(dates[previ]):
            previ = -1

        # End of an interval -> save
        if previ >= 0 and date == t_iter.end(date):
            tslc.append(slice(previ,i+1))
            cnt_out.append(i+1 - previ)
            start_dates.append(dates[previ])
            if i+1 >= len(dates):
                end_dates.append(dates[-1] + dtd)
            else:
                end_dates.append(dates[i+1])
            previ = -1
        
        # Start of an interval -> remember 
        if date == t_iter.start(date):
            previ = i
    
    outlen = len(tslc)

    # 2. Initialising the output array
    shape = list(dat.shape)
    shape[0] = outlen
    shape = tuple(shape)

    dat_out = np.empty(shape)
    if not type(bins) == type(None):
        valid_out = np.empty(shape[0:1] + (len(bins),) + shape[1:], dtype='i4')
    else:
        valid_out = np.empty(shape, dtype='i4')
    
    # 3. Doing the actual calculations
    for i in range(outlen):
        dat_ = dat[tslc[i]]
        if not type(bins) == type(None):
            for bi in range(len(bins)-1):
                upper, lower = bins[bi+1], bins[bi]
                if upper > lower:
                    valid_out[i,bi,:,:] += np.logical_and(dat_ >= lower, dat_ <  upper).sum(axis=0)
                else:
                    valid_out[i,bi,:,:] += np.logical_or(dat_ <  upper, dat_ >= lower).sum(axis=0)
            dat_out[i] = cal_mfv(valid_out[i,:,:,:], bins)

        else:
            dat_out[i] = np.nanmean(dat_, axis=0)
            valid_out[i] = cnt_out[i] - np.isnan(dat_).sum(axis=0)

    return start_dates, end_dates, dat_out, valid_out, np.array(cnt_out)



#
# Varimax rotation for EOFs as introduced in Kaiser (1958)
#  -> if raw = False use normal varimax, otherwise raw varimax
def varimax(phi, raw=False, gamma = 1.0, max_iter = 100, tol = 1e-5):
    ''' Varimax rotation for EOFs

    As introduced by Kaiser (1958).

    Parameters
    ----------
    phi : np.ndarray
        EOF loadings
    raw : bool
        Optional, default ``False``. If ``True`` use the "raw varimax" rotation instead of the "normal varimax".
    gamma : float
        Optional, default ``1``. Iteration parameter.
    max_iter : int
        Optional, default ``20``. Maximum number of iterations to find the optimal rotation.
    tol : float
        Optional, default ``1.0e-6``. Numerical tolerance to consider the iteration converged.
        
    Returns
    -------
    np.ndarray
        Rotated EOF loadings
    np.ndarray
        Rotation matrix
    '''

    if not raw:
        norm = np.sqrt((phi**2).sum(axis=1))
        phi /= norm[:,np.newaxis]
    
    p, k = phi.shape
    R = np.eye(k)
    d = 0
    for i in range(max_iter):
        Lambda = np.dot(phi, R)
        u,s,vh = sp.linalg.svd(np.dot(
            phi.T, 
            np.asarray(Lambda)**3 - (gamma/p) * np.dot(Lambda, np.diag(np.diag(np.dot(Lambda.T,Lambda))))
        ))
        R = np.dot(u,vh)
        d_old = d
        d = np.sum(s)
        if d < d_old * (1 + tol):
            break
        if i == max_iter - 1:
            print('WARNING: Max iteration reached without convergence; try fiddeling with tolerance and maximum number of iterations.')
    
    rphi = np.dot(phi, R)

    if not raw:
        rphi *= norm[:,np.newaxis]

    return rphi, R



def dist_sphere(lon1, lat1, lon2, lat2, r=6.37e6):
    ''' Shortest distance on a sphere

    Calculate the great circle distance between two points on the surface of a sphere, 
    using spherical trigonometry. By default, the radius of the sphere is assumed to 
    be the Earth radius, R = 6370 km, but that can be changed via the optional 
    parameter r.

    Both the first and second points can be an array of points. If both points are
    actually arrays of points, these arrays must have compatible shapes in the sense of 
    the numpy broadcasting rules.

    Parameters
    ----------
    lon1 : float or np.ndarray
        Longitude(s) of the first point(s) in degrees.
    lat1 : float or np.ndarray
        Latitude(s) of the first point(s) in degrees.
    lon2 : float or np.ndarray
        Longitude(s) of the second point(s) in degrees.
    lat2 : float or np.ndarray
        Latitude(s) of the second point(s) in degrees.
    r : float or np.ndarray
        *Optional*. Radius of the sphere(s). Defaults to the Earth radius.
    
    Returns
    -------
    float or np.ndarray
        Distance(s) between the first and second points
    '''
        
    dlon = np.pi/180 * (lon2 - lon1)
    lat1r = np.pi/180 * lat1
    lat2r = np.pi/180 * lat2
    acos = np.sin(lat1r)*np.sin(lat2r) + np.cos(lat1r)*np.cos(lat2r)*np.cos(dlon)
    dist = r * np.arccos(np.maximum(np.minimum(acos,1.0),-1.0))

    return dist



def dist_green_latlon(lon, lat):
    ''' Calculate distances from points along zero meridian and 0°-90°N

    These distances can be used to readily look up the distance between any two points 
    on the given latitude/longitude raster.

    This approach is not unlike that of Green's functions to solve PDEs, hence the name.

    Parameters
    ----------
    lon : list or np.ndarray with dimensions (nx)
        Longitudes of grid points
    lat : list or np.ndarray with dimensions (ny)
        Latitudes of grid points
    
    Returns
    -------
    np.ndarray with dimensions (nx//2,ny,nx)
        Distances between all points and points along the zero meridian between 0° and 90°N in meters
    '''

    if not 0 in lon:
        raise ValueError('Zero meridian must be part of the given longitudes')

    Ngreen = len(lat)//2 + 1
    dists = np.empty((Ngreen, len(lat), len(lon)))
    
    lon = np.array(lon)
    lat = np.array(lat)

    for n in range(Ngreen):
        dists[n,:,:] = dist_sphere(lon[np.newaxis,:], lat[:,np.newaxis], 0, lat[n])
    
    return dists



def dist_from_mask_latlon(featmask, green_dists, izero, quiet=True):
    ''' Calculate the minimum distance to any feature, as given by a mask field

    Parameters
    ----------
    featmask : np.ndarray with dimensions (...,ny,nx) and dtype bool
        Binary mask field marking the locations (by value True) of detected features.
    green_dists : np.ndarray with dimensions (nx//2+1,ny,nx)
        Precalculated array of distances between any two points on the grid. Use 
        ``dist_green_latlon`` to prepare this array.
    izero : int
        Longitude-index of the zero meridian
    quiet : bool
        *Optional*, default ``True``. If not quiet, progress is printed.
    
    Returns
    -------
    np.ndarray with dimensions (...,ny,nx)
        Minimum distance to any detected feature in meters.
    '''

    ny = featmask.shape[-2]
    mindists = np.ones(featmask.shape) * 50.0e6 # Larger than the circumpherence of the Earth
    poss = np.argwhere(featmask)
    updateevery = len(poss) // 1000
    for idx, pos in enumerate(poss):
        if idx % updateevery == 0 and not quiet:
            print(f'{idx/10/updateevery:3.1f}%', end=chr(13))

        overwrite = tuple(pos[:-2]) + (slice(None), slice(None))
        j, i = pos[-2:]
        shift = i - izero
        if j < green_dists.shape[0]:
            curdists = np.roll(green_dists[j,:,:], shift, axis=1)
        else:
            j = ny - j - 1
            curdists = np.roll(green_dists[j,::-1,:], shift, axis=1)

        mindists[overwrite] = np.minimum(mindists[overwrite], curdists)

    if not quiet:
        print('')
    
    return mindists



def smooth_xy_nan(dat, nsmooth):
    ''' As smooth_xy, but first using the fill_nan function to remove potential NaNs.

    The NaN mask is conserved and applied to the smoothed field before it is returned.

    Parameters
    ----------
    dat : np.ndarray with dimensions (nz,ny,nx)
        Input data to be smoothed.
    nsmooth : int
        Number of passes of the three-point filter
    
    Returns
    -------
    np.ndarray with dimensions (nz,ny,nx)
        Smoothed version of the input data.
    '''
    
    # Save the NaN mask
    mask = np.isnan(dat)
    
    # Remove zonal-time mean before filling the NaNs to make the initial guess of zeros a bit more applicable
    zonaltimemean = np.nanmean(dat, axis=2).mean(axis=0)
    dat_ = dat - zonaltimemean[np.newaxis,:,np.newaxis]
    dat_, conv = fill_nan(dat_)
    dat_ += zonaltimemean[np.newaxis,:,np.newaxis]
    
    # Smooth the filled field and then apply the NaN mask again
    dats = smooth_xy(dat_, nsmooth)
    dats[mask] = np.nan

    return dats



def lanczos_weights_lowpass(cutoff, window=None, size=2):
    ''' Calculate the weights for a 1d Lanczos lowpass filter with a given cutoff

    The function is following Duchon C. E. (1979): Lanczos Filtering in One and Two 
    Dimensions. Journal of Applied Meteorology, Vol 18, pp 1016-1022. It is adapted 
    from code found at several places on the internet, so proper attribution is 
    unclear. Sources were

     - https://scitools.org.uk/iris/docs/v1.2/examples/graphics/SOI_filtering.html
     - https://github.com/liv0505/Lanczos-Filter

    but the code might live in and originate from further unidentified places.

    The window can be either specified manually, or can be calculated to achieve
    a filter of defined size (often called "a"). If nothing else is specified a 
    size-2 filter is constructed (one positive central lobe with two negative lobes 
    on either side).

    Parameters
    ----------
    cutoff : int/float
        The cutoff frequency in time steps.
    window : int
        *Optional*, by default calculated to yield a filter of given order. The length 
        of the filter window in time steps.
    size : int
        *Optional*, default 2. If the window length is calculated automatically, the
        size parameter of the requested filter can be specified here. The resulting 
        window length is ``int(size*cutoff) + 1``.

    Returns
    -------
    np.ndarray with dimension (window)
        The filter weights.
    '''

    if not window:
        window = int(size*cutoff) + 1
    
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 / cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi / cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma

    # Normalize to sum of 1
    w[1:-1] /= w[1:-1].sum()
    
    return w[1:-1]

#
