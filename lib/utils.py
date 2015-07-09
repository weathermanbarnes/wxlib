#!/usr/bin/env python

''' Utilities that don't fit elsewhere '''

import sys

import dynfor
import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.utils, 'utils', sys.modules[__name__])


import math
import numpy as np
from scipy.special import erfinv

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
	    With the default data layout, only relevant data needs to be read when only 
	    a time slice of the entire data is requested. Hence, using cut to limit your 
	    data request can make reading the data largely more efficient.
	
	Returns
	-------
	np.ndarray
	    The scaled data contained in the nc.NetCDFVariable object `var`.
	"""
	
	# Apply cut first, for speed
	dat = var[cut,::]
	dat = dat.astype('f8')

	# Mask missing and fill values
	if hasattr(var, 'missing_value'):
		dat[dat == var.missing_value] = np.nan
	if hasattr(var, '_FillValue'):
		dat[dat == var._FillValue] = np.nan
	
	# Apply scaling, numpy is faster than Fortran function!
	if hasattr(var, 'scale_factor') or hasattr(var, 'add_offset'):
		dat = dat*var.scale_factor + var.add_offset
	
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
	
	missing = -32767

	maxv  = np.nanmax(var)
	minv  = np.nanmin(var)

	# divide in 2^16-2 intervals, values from -32766 -> 32767 ; reserve -32767 as missing value
	scale = (maxv-minv)/65534.0
	off   = +32766.5*scale + minv

	res = np.round((var[::] - off)/scale)

	res[np.isnan(res)] = missing

	return res.astype('>i2'), scale, off, missing


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
		raise NotImplementedError, 'Concatenation not implemented for %d dimensions' % len(data.shape)
	
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

#
# return a 3d boolean array where frontal points are True, elsewere False
def mask_fronts(fronts, froff, s=(361,720)):
	''' To be made obsolete by saving cold/warm/stat fronts separately as lines in the standard-dynlib way 
	
	See also
	--------
	``mask_lines``, ``smear_lines``.
	'''
	masks = [np.zeros((len(fronts), s[0], s[1]), dtype='bool'),
		 np.zeros((len(fronts), s[0], s[1]), dtype='bool'),
		 np.zeros((len(fronts), s[0], s[1]), dtype='bool') ]

	for t in range(len(fronts)):
		for f in range(3):
			for n in range(froff[t,f].max()):
				# python starts counting at zero, unlike fortran
				j = round(fronts[t,f,n,1] -1)
				i = round(fronts[t,f,n,0] -1) % s[1]
				masks[f][t,j,i] = True

	return masks


def mask_lines_with_data(lines, loff, dat=None, s=(361,720)):
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
	dat : np.ndarray with dimensions (y,x)
	    Optional: Data to be used for marking on the map
	s : 2-tuple of int
	    Optional: Grid dimensions
	
	Returns
	-------
	np.ndarray
	    Gridded map of lines
	'''

	mask = np.zeros((lines.shape[0], s[0], s[1]))

	for t in range(lines.shape[0]):
		for n in range(loff[t].max()):
			# python starts counting at zero, unlike fortran
			j = round(lines[t,n,1] -1)
			i = round(lines[t,n,0] -1) % s[1]
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


def smear_lines(lines, loff, s=(361,720)):
	''' Mask lines in a gridded map and then smooth slightly

	Grid points containing a line will be marked with the value ``1``,
	otherwise the retuned map contains zeros.

	Notes
	-----
	
	 * Grid size should become a parameter in the configuration
	 * Make the filter configurable

	Parameters
	----------
	lines : np.narray with dimensions (pointindex,infotype)
	    Lines to be marked on the map
	loff : np.array with dimensions (lineindex)
	    List of point indexes for the first points of each line
	s : 2-tuple of int
	    Optional: Grid dimensions
	
	Returns
	-------
	np.ndarray
	    Gridded map of lines
	'''
	
	filtr_len = 5
	filtr_func = mk_gauss(0, 1)
	filtr = np.array(map(filtr_func, range(-filtr_len,filtr_len+1)))
	filtr /= sum(filtr)

	mask = dynfor.utils.mask_lines(s[1], s[0], lines, loff)
		
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

	retlon = []
	retlat = []
	retxy  = []
	prevxy = 0
	for i, (lon1, lat1) in zip(range(len(coords)-1), coords[:-1]):
		lon2, lat2 = coords[i+1]

		(x1,x2), (y1,y2) = m((lon1,lon2), (lat1,lat2))
		d12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
		N = np.ceil(d12/dxy)

		x  = [x2*alpha + x1*(1-alpha) for alpha in np.arange(0.0,(N+1)/N,1.0/N)]
		y  = [y2*alpha + y1*(1-alpha) for alpha in np.arange(0.0,(N+1)/N,1.0/N)]
		xy = [d12*alpha + prevxy for alpha in np.arange(0.0,(N+1)/N,1.0/N)]

		prevxy = xy[-1]

		lon, lat = m(x, y, inverse=True)
		retlon.extend(lon)
		retlat.extend(lat)
		retxy.extend(xy)

	return retlon, retlat, retxy

def aggregate(dates, dat, agg):
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
	    * ``met_season``: Seasons after their standard meteorological definition
	    * ``cal_monthly``: Months following the calendar definition
	    * ``10d``: 10-day intervals starting from the first given date
	    * ``cal_weekly``: Weeks following the calendar definition
	    * ``7d`` or ``weekly``: 7-day intervals starting from the first given date, 
	      irrespective of the day of the week
	    * ``pendad``: Pentads according to their definition
	    * ``5d``: 5-day intervals starting from the first given date
	    * ``3d``: 3-day intervals starting from the first given date
	    * ``2d``: 2-day intervals starting from the first given date
	    * ``1d`` or ``daily``: 1-day intervals starting from the first given date
	
	Returns
	-------
	list of datetime
	    Dates at with the aggregation periods start
	np.ndarray
	    Aggregated data
	'''

	dtd = dates[1] - dates[0]

	# Helper functions generating functions :-)
	def first_func_gen(tdays):
		epoch = dates[0]
		return lambda date: epoch + td(((date - epoch).days / tdays)*tdays)
	def last_func_gen(tdays):
		epoch = dates[0]
		return lambda date: epoch + td(((date - epoch).days / tdays + 1)*tdays) - dtd

	# 1a. Defining aggregation types 
	tslc = []
	dates_out = []
	
	# Aggregate everything
	if agg == 'all':
		first_func = lambda date: dates[0]
		last_func = lambda date: dates[-1]
		agg = 'func'
	
	# Aggregate by season
	elif agg == 'met_season':
		# Defining some helpers
		season_start = {1: 12, 2: 12, 3: 3, 4: 3, 5: 3, 6: 6, 7: 6, 8: 6, 9: 9, 10: 9, 11: 9, 12: 12}
		first_year_offset = {1: -1, 2: -1, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0}
		last_year_offset = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 1}

		first_func = lambda date: dt(date.year + first_year_offset[date.month], season_start[date.month], 1)
		last_func = lambda date: dt(date.year + last_year_offset[date.month], (season_start[date.month]+2) % 12 + 1, 1) - dtd
		agg = 'func'
	
	# Aggregate by calendar months
	elif agg == 'cal_monthly':
		first_func = lambda date: dt(date.year,date.month,1)
		last_func = lambda date: dt(date.year+date.month/12, date.month % 12 + 1, 1) - dtd
		agg = 'func'
	
	# Aggregate by ISO calendar weeks
	elif agg == 'cal_weekly':
		# date.isocalendar gives the a tuple (ISOWEEK_YEAR, ISOWEEK, ISOWEEK_DAY)
		first_func = lambda date: date - (date.isocalendar()[2] - 1)*td(1,0) - td(0,3600*date.hour)
		last_func = lambda date: date + (8 - date.isocalendar()[2])*td(1,0) - td(0,3600*date.hour) - dtd
		agg = 'func'
	
	# Aggegate by pendats, following the CPC definition that Pentad 12 is always 25 Feb -- 1 Mar, even for leap years
	elif agg == 'pentad':
		def pentad_off(date):
			leap = 0
			startleap = 0
			if calendar.isleap(date.year) and date.timetuple().tm_yday >= 60:
				leap = 1
			pentad = ((date - dt(date.year, 1, 1)).days - leap)/ 5 
			if leap > 0 and pentad > 11:
				startleap = 1
			return pentad, startleap, leap

		def first_func(date):
			pentad, startleap, endleap = pentad_off(date)
			return dt(date.year, 1, 1) + pentad*td(5) + td(startleap)
		def last_func(date):
			pentad, startleap, endleap = pentad_off(date)
			return dt(date.year, 1, 1) + (pentad+1)*td(5) + td(endleap) - dtd
		agg = 'func'
	
	# Aggregate by 10-day periods
	elif agg == '10d':
		first_func = first_func_gen(10)
		last_func = last_func_gen(10)
		agg = 'func'
	# Aggregate by 7-day periods
	elif agg == '7d' or agg == 'weekly':
		first_func = first_func_gen(7)
		last_func = last_func_gen(7)
		agg = 'func'
	# Aggregate by 5-day periods
	elif agg == '5d':
		first_func = first_func_gen(5)
		last_func = last_func_gen(5)
		agg = 'func'
	# Aggregate by 3-day periods
	elif agg == '3d':
		first_func = first_func_gen(3)
		last_func = last_func_gen(3)
		agg = 'func'
	# Aggregate by 2-day periods
	elif agg == '2d':
		first_func = first_func_gen(2)
		last_func = last_func_gen(2)
		agg = 'func'
	# Aggregate by 24h periods
	elif agg == '1d' or agg == 'daily':
		first_func = first_func_gen(1)
		last_func = last_func_gen(1)
		agg = 'func'
	
	# Unknown aggregation period
	else:
		raise NotImplementedError, 'Unknown aggregation specifier `%s`' % str(agg)
	
	# 1b. Estimating length of output array
	#     Aggegate by output of functions first_func and last_func; 
	#     These functions give the first and the last time step for the interval `date` belongs to
	if agg == 'func':
		previ = -1
		for i, date in zip(range(len(dates)), dates):
			# Previous interval unfinished, yet we are in a new one.
			# -> We apparently jumped over a range of dates and have to find a new start
			if previ >= 0 and not first_func(date) == first_func(dates[previ]):
				previ = -1
			
			# End of an interval -> save
			if previ >= 0 and date == last_func(date):
				tslc.append(slice(previ,i))
				dates_out.append(dates[previ])
				previ = -1
			
			# Start of an intercal -> remember 
			if date == first_func(date):
				previ = i
	
	else:
		raise ValueError, 'Unknown agg value `%s`' % str(agg)
	
	outlen = len(tslc)

	# 2. Initialising the output array
	s = list(dat.shape)
	s[0] = outlen
	s = tuple(s)

	dat_out = np.empty(s)
	
	# 3. Doing the actual calculations
	for i in xrange(outlen):
		dat_out[i] = dat[tslc[i]].mean(axis=0)

	return dates_out, dat_out



#
# Varimax rotation for EOFs as introduced in Kaiser (1958)
#  -> if raw = False use normal varimax, otherwise raw varimax
def varimax(phi, raw=False, gamma = 1.0, max_iter = 20, tol = 1e-6):
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
	for i in xrange(max_iter):
		d_old = d
		Lambda = np.dot(phi, R)
		u,s,vh = np.linalg.svd(np.dot(
			phi.T, 
			np.asarray(Lambda)**3 - (gamma/p) * np.dot(Lambda, np.diag(np.diag(np.dot(Lambda.T,Lambda))))
		))
		R = np.dot(u,vh)
		d = np.sum(s)
		if d/d_old < tol: break
	
	rphi = np.dot(phi, R)

	if not raw:
		rphi *= norm[:,np.newaxis]

	return rphi, R

#
