#!/usr/bin/env python
# -*- encoding: utf-8

''' Input/Output library for meteorological data

Provides convenient access to met-data stored in numpy, netCDF and Matlab formats,
and provides functions to save meteorological fields and lines in netCDF files.
'''

import os
import sys

import numpy as np
import netCDF4 as nc
import scipy.io.matlab as mat
import pytz
import calendar

from settings import conf
import utils
from gridlib import grid_by_nc, grid_by_static

from datetime import datetime as dt, timedelta as td
import tagg

import dynfor
dynlib_version = (''.join(dynfor.consts.version)).strip()




# #############################################################################
# 1. Open all sorts of files
# 

# TODO: Restrict the expected file types by providing a file name extension
def metopen(filename, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False, no_static=False):
	''' Find and open files by name
	
	Uses the conf.datapath list to locale files in a variety of different locations.
	The located files might either be numpy-files, netCDF-files or matlab mat-files. 
	For each of the file types, metopen returns the requested variable and some meta-
	information about the variable, if not suppressed.

	Parameters
	----------
	filename : str
	    The name of the file, excluding the file ending.
	q : str
	    *Optional*. The requested variable within the file.
	cut : slice
	    *Optional*, default ``slice(None)``. Limit the request to a given time slice. 
	    With the default data layout, only relevant data needs to be read when only 
	    a time slice of the entire data is requested. Hence, using cut to limit your 
	    data request can make reading the data largely more efficient.
	verbose : bool
	    *Optional*, default ``False``. Print debug information on which files are 
	    being looked for.
	no_dtype_conversion : bool
	    *Optional*, default ``False``. By default, ``metopen`` uncompresses data in the 
	    file and converts all data to float64. This behaviour can be suppressed by 
	    setting this option to ``True``. This implies, however, that scaling and offset
	    cannot be applied automatically.
	no_static : 
	    *Optional*, default ``False``. By default, ``metopen`` does its best to 
	    provide meta-information about the requested file, using 
	    :module:`grid.gridlib`, and returns the meta-information as a third value. 
	    This behaviour can be suppressed by setting this parameter to ``True``.
	
	Returns
	-------
	data file object
		python data, netCDF or Matlab file object.
	np.ndarray
		If q given, data of the requested variable.
	grid.gridlib
		If ``no_static=False`` meta-information about the requested data.
	'''
	
	if not type(cut) == slice:
		raise ValueError, 'cut must be a 1-dimensional slice object' 
	
	tried = []
	for path in conf.datapath:
		static = None

		if verbose:
			print 'Trying: '+path+'/'+filename+'.*'

		if os.path.exists(path+'/'+filename+'.npy'):
			if q:
				dat = np.load(path+'/'+filename+'.npy', mmap_mode='r')
				dat = dat[cut]
			print 'Found '+path+'/'+filename+'.npy'
			f = None
		elif os.path.exists(path+'/'+filename+'.npz'):
			f = np.load(path+'/'+filename+'.npz')
			if q:
				if q not in f.files:
					tried.append(path+'/'+filename+'.npz')
					continue
				dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.npz'
		elif os.path.exists(path+'/'+filename+'.mat'):
			f = mat.loadmat(path+'/'+filename+'.mat')
			if q:
				if q not in f:
					tried.append(path+'/'+filename+'.mat')
					continue
				dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.mat'
		elif os.path.exists(path+'/'+filename+'.nc'):
			f   = nc.Dataset(path+'/'+filename+'.nc', 'r')
			if q:
				var = f.variables[q]
				if q not in f.variables:
					tried.append(path+'/'+filename+'.nc')
					continue
				if not no_dtype_conversion:
					dat = utils.scale(var, cut=cut)
				else:
					dat = var[cut]
			else:
				var = None

			if not no_static:
				static = grid_by_nc(f, var)
				# TODO: Where to search for topography in nc files?
				static.oro = np.zeros((static.ny, static.nx))

			print 'Found '+path+'/'+filename+'.nc'
		else:
			tried.append(path)
			continue
		
		if q and not no_dtype_conversion:
			dat = dat.astype('f8')

		if not no_static:
			if not static:
				static = get_static(verbose, no_dtype_conversion)
		else:
			if q:
				return f, dat
			else:
				return f
		if q:
			return f, dat, static
		else:
			return f, static
	
	raise ValueError, '%s.* not found in any data location. \nTried the following (in order):\n\t%s' % (filename, '\n\t'.join(tried))


# Save dat as a netCDF file, using the metadata in static. 
def metsave(dat, static, q, plev, compress_to_short=True):
	''' Save data in a netCDF file

	Parameters
	----------
	dat : np.ndarray
		Data to be saved.
	static : gridlib.grid
		Some meta information about the data, like the grid information.
	q : str
		The variable name identifier, following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``.
	plev : str
		String representation of the vertical level on which the data is defined, e.g. ``'700'`` 
		for 700 hPa or ``'pv2000'`` for the PV2-surface.
	compress_to_short : bool
		*Optional*, default ``True``. By default, ``metsave`` compresses the data by converting
		the data field into int16, using the float64 ``add_offset`` and ``scale_factor`` attributes 
		to represent the data.
	'''
	s = dat.shape
	if not len(s) == 3 or (conf.gridsize and not s[1:] == conf.gridsize):
		raise NotImplementedError, 'dat does not seem to be a context-conform (t,y,x)-array.'
	
	if not conf.oformat == 'nc':
		raise NotImplementedError, 'Currently only saving in netCDF implemented in metsave.'

	if not s[0] == len(static.t):
		raise ValueError, 'Time dimension in data (%s) and static (%s) are not equally long.' % (s[0], len(static.t))
	if not s[1] == len(static.y[:,0]):
		raise ValueError, 'y-dimension in data (%s) and static (%s) are not equally long.' % (s[1], len(static.y[:,0]))
	if not s[2] == len(static.x[0,:]):
		raise ValueError, 'x-dimension in data (%s) and static (%s) are not equally long.' % (s[2], len(static.x[0,:]))

	now = dt.now(pytz.timezone(conf.local_timezone))
	of = nc.Dataset(conf.opath+'/'+(conf.file_std % {
		'time': dts2str(static.t_parsed), 'plev': plev, 'qf': conf.qf[q]})+'.nc', 'w', format='NETCDF3_CLASSIC')
	of.setncatts({'Conventions': 'CF-1.0', 
			'history': '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	})

	of.createDimension('time', s[0])
	of.createDimension(static.x_name, s[2])
	of.createDimension(static.y_name, s[1])

	ot = of.createVariable('time', 'i', ('time',))
	ot.setncatts({'long_name': 'time', 'units': static.t_unit})
	ot[::] = static.t
	olat = of.createVariable(static.y_name, 'f', (static.y_name,))
	olat.setncatts({'long_name': static.y_name, 'units': static.y_unit})
	olat[::] = static.y[:,0]
	olon = of.createVariable(static.x_name, 'f', (static.x_name,))
	olon.setncatts({'long_name': static.x_name, 'units': static.x_unit})
	olon[::] = static.x[0,:]
	
	if compress_to_short:
		ovar = of.createVariable(q, 'h', ('time', static.y_name, static.x_name,))
		dat, scale, off, fill = utils.unscale(dat)
		ovar.setncatts({'long_name': conf.q_long[q], 'units': conf.q_units[q],
				'add_offset': off, 'scale_factor': scale})
		if fill: 
			ovar._FillValue = fill
			ovar.missing_value = fill
	else:
		ovar = of.createVariable(q, 'f', ('time', static.y_name, static.x_name,))
		ovar.setncatts({'long_name': conf.q_long[q], 'units': conf.q_units[q]})
	ovar[::] = dat

	of.close()

	return


# Save dat as a netCDF file, using the metadata in static. 
def metsave_lines(dat, datoff, static, time, plev, q, qoff):
	''' Save line data in a netCDF file

	A disk-space-effective way to store collections of lines (e.g. front lines, jet axes, ...)
	in a netCDF file. The coordinates of all points belonging to a line are stored in a long
	list in ``dat``. The beginnings and endings of lines are defined by ``datoff``, which lists
	the starting point indexes of all lines in the collection. 

	Parameters
	----------
	dat : np.ndarray with dimensions (time,pointindex,infotype)
		Coordinates and an additional piece of information about each point belonging to 
		a line.
	datoff : np.ndarray with dimensions (time,lineindex)
		Starting point indexes for each line.
	static : gridlib.grid
		Some meta information about the data, like the grid information.
	time : str
		String representation of the time period covered by the data, e.g. ``'2005'`` for the year
		2005 or ``'20050317'`` for 17 March 2005.
	plev : str
		String representation of the vertical level on which the data is defined, e.g. ``'700'`` 
		for 700 hPa or ``'pv2000'`` for the PV2-surface.
	q : str
		The variable name identifier, e.g. ``'jetaxis'``.
	qoff : str
		The variable name identifier for the list of offsets, e.g. ``'jaoff'``. 
	'''
	if not len(dat.shape) == 3 and len(datoff.shape) == 2:
		raise RuntimeError, 'dat and/or datoff have the wrong number of dimensions'
	if not dat.shape[0] == datoff.shape[0]:
		raise RuntimeError, 'dat and datoff have different time axes'
	if not dat.shape[2] == 3:
		raise RuntimeError, 'dat does not have size 3 in the third dimension'

	now = dt.now(pytz.timezone('Europe/Oslo'))
	of = nc.Dataset(conf.opath+'/'+(conf.file_std % {'time': time, 'plev': plev, 'qf': conf.qf[q]})+'.nc', 'w', format='NETCDF3_CLASSIC')
	of.setncatts({'Conventions': 'CF-1.0', 
			'history': '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	})
	
	# The maximum amount of line points for all time steps
	#llen = int(datoff.max())
	llen = dat.shape[1] # Having diffent dimensions for each file makes reading much harder
	# The maxmimum amount of lines for all time steps
	#olen = int(np.max([datoff_t.argmin() for datoff_t in datoff[:,1:]]))
	olen = datoff.shape[1]

	of.createDimension(static.t_name, dat.shape[0])
	of.createDimension('pointindex', llen)
	of.createDimension('infotype', 3)
	of.createDimension('lineindex', olen)
	
	# Define helper dimensions; not explicitly used in any variable but useful to interpret the grid point indexes stored in dat
	of.createDimension(static.x_name, static.x.shape[1])
	of.createDimension(static.y_name, static.y.shape[0])
	olon = of.createVariable(static.x_name, 'f', (static.x_name,))
	olon.setncatts({'long_name': static.x_name, 'units': static.x_unit})
	olon[:] = static.x[0,:]
	olat = of.createVariable(static.y_name, 'f', (static.y_name,))
	olat.setncatts({'long_name': static.y_name, 'units': static.y_unit})
	olat[:] = static.y[:,0]

	ot = of.createVariable(static.t_name, 'i', ('time',))
	ot.setncatts({'long_name': static.t_name, 'units': static.t_unit})
	ot[:] = static.t
	olidx = of.createVariable('pointindex', 'i', ('pointindex',))
	olidx.setncatts({'long_name': 'Index of point along all lines', 'units': '1'})
	olidx[:] = range(llen)
	olity = of.createVariable('infotype', 'c', ('infotype',))
	olity.setncatts({'long_name': 'Type of info stored for point', 'units': 'enum'})
	olity[:] = ['X', 'Y', 'I']
	ooidx = of.createVariable('lineindex', 'i', ('lineindex',))
	ooidx.setncatts({'long_name': 'Index of line', 'units': '1'})
	ooidx[:] = range(olen)
	
	oq = of.createVariable(q, 'f', ('time', 'pointindex', 'infotype',))
	oq.setncatts({'long_name': conf.q_long[q], 'units': 'mixed'})
	oq[::] = dat[:,:llen,:]

	oqoff = of.createVariable(qoff, 'i', ('time', 'lineindex',))
	oqoff.setncatts({'long_name': 'Index of first point of line', 'units': '1'})
	oqoff[::] = datoff[:,:olen]

	of.close()

	return


# Get static information
def get_static(verbose=False, no_dtype_conversion=False):
	''' Get standard meta-information (for ERA-Interim)

	Parameters
	----------
	verbose : bool
		*Optional*, default ``False``. Print debug information on where the static
		file is being looked for.
	no_dtype_conversion : bool
		*Optional*, default ``False``. By default, ``metopen`` uncompresses data in the 
		file and converts all data to float64. This behaviour can be suppressed by 
		setting this option to ``True``.

	Returns
	-------
	'''
	fo, oro = metopen(conf.file_static, 'oro', verbose=verbose, no_dtype_conversion=no_dtype_conversion, no_static=True)
	static = grid_by_static(fo)
	static.oro = oro[::]
	fo.close()

	return static


# #############################################################################
# 2. Generalised data fetchers
# 

def get_instantaneous(q, dates, plevs=None, tavg=False, force=False, **kwargs):
	''' Data fetcher for instantaneous or short-term averaged fields

	Allows general data requests in the configured data base, e.g. ERA-Interim. The request
	can span several files, e.g. by including several vertical levels or by covering several
	years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
	However, length-1 dimensions are removed, so the returned array can have fewer 
	dimensions.

	Parameters
	----------
	q : str
		Variable name identifier, following the ECMWF conventions as far as appicable, 
		e.g. ``'u'`` or ``'msl'``.
	dates : list of datetime
		The minimum and maxmimum dates in this list define the requested time interval.
	plevs : list of str or str
		*Optional*, defaults to all available pressure levels. The string representation
		of the requested vertical level(s), e.g. ``'700'`` for 700 hPa or ``'pv2000'`` 
		for the PV2-surface.
	tavg : bool
		*Optional*, default ``False``. Instead of returning a time-dimension, return
		the temporal average of the requested data.
	force : bool
		*Optional*, default ``False``. Turn off the error, if large amounts of data are
		requested at once. **Be sure you know what you are doing, when setting this to 
		``True``! Your request might make your script occupy a large fraction of the 
		system memory**.
	
	Keyword arguments
	-----------------
	metopen arguments : all optional
		Optional arguments passed on to calls of metopen within this function.
	'''
	
	if not plevs:
		plevs = conf.plevs
	elif not type(plevs) == np.ndarray and not type(plevs) == list:
		plevs = [plevs,]
	
	if conf.years:
		dt0 = dt(conf.years[0],1,1,0)
	elif conf.times: 
		dt0 = conf.times[0]
	else:
		raise ValueError, 'Could not determine start date of the data set, either conf.years or conf.times must be provided.'

	dts = conf.timestep.total_seconds()

	# Convert dates to time indexes
	if type(dates) not in ([np.ndarray, list, tuple, set]):
		dates = [dates, ]
	years = set(map(lambda date: date.year, dates))
	years = np.arange(min(years), max(years)+1)
	tidxs = map(lambda date: date-dt0, dates)
	tidxs = map(lambda diff: int(diff.total_seconds()/dts), tidxs)
	tsmin, tsmax = min(tidxs), max(tidxs)
	
	# Checking max length
	maxtlen = 3000
	if tsmax-tsmin > maxtlen:
		if force:
			print 'Warning: you requested %d time steps.' % (tsmax-tsmin)
		else:
			raise ValueError, 'Cowardly refusing to fetch %d > %d time steps.\nUse force=True to override.' % (tsmax-tsmin, maxtlen)
	
	# Remove no_static if present
	kwargs.pop('no_static', None)
	static = None

	dat = None
	for year in years:
		# Construct the slice
		fst = dt(year,1,1,0) - dt0
		lst = (dt(year+1,1,1,0) - conf.timestep) - dt0
		fst = int(fst.total_seconds()/dts)
		lst = int(lst.total_seconds()/dts) + 1

		# Leave out unnecessary indexes for better compatibility
		cut = slice(max(tsmin - fst, 0),min(1+tsmax - fst, lst - fst))
		datcut = slice(fst+cut.start-tsmin, fst+cut.stop-tsmin)
		
		# One or more vertical levels?
		i = 0
		for plev in plevs:
			if type(dat) == type(None):
				f, d, static = metopen(conf.file_std % {'time': year, 'plev': plev, 'qf': conf.qf[q]}, q, cut=cut, **kwargs)
				s = (1+tsmax-tsmin, len(plevs), ) + d.shape[1:]
				dat = np.empty(s, dtype=d.dtype)
				static.t = static.t[cut]
				if type(static.t_parsed) == np.ndarray:
					static.t_parsed = static.t_parsed[cut]
			else:
				f, d, static_ = metopen(conf.file_std % {'time': year, 'plev': plev, 'qf': conf.qf[q]}, q, cut=cut, **kwargs)
				static.t = np.concatenate((static.t, static_.t[cut]))
				if type(static.t_parsed) == np.ndarray:
					static.t_parsed = np.concatenate((static.t_parsed, static_.t_parsed[cut]))
			
			dat[datcut,i,::] = d
			i += 1

	# Time-averaging if specified
	if tavg and len(dates) > 1:
		dat = dat.mean(axis=0)
	
	dat = dat.squeeze()
	
	return dat, static


# Get aggregated (average, standard deviation, etc.) fields
def get_aggregate(q, dates, plevs=None, tavg=False, force=False, **kwargs):
	''' Request aggregate data 
	
	This function is not implemented yet. It is only included to define the future API.
	'''
	raise NotImplementedError, 'If you need it, implement it!'

	return dat


# #############################################################################
# 3. Utilities
#

def dts2str(dates, agg=None):
	''' Find a short, but descriptive string representation for a list of dates

	If only a single date is given, the formatted date of the form YYYYMMDDHH is returned.

	If several dates are given, the minimum and maximum dates are returned. The format 
	depends on the dates. For the start date, hours, days and month are potentially
	omitted if they are the first in their respective superordinate period. Analogously, 
	hours, days and month are omitted for the last date if they correspond to the last 
	time step within their superordinate period. 
	
	The calculation of the last time step depends on the ``conf.timestep`` property.
	For example the last time step of October 2014 in the ERA-Interim data is 18 UTC
	on 31 October 2014.

	If the given number of dates given matches the number of dates expected for the given 
	start and end dates with the ``conf.timestep``, the dates are *assumed* to be
	contiguous, and the two dates are joined by a minus sign ``'-'``. Otherwise, the 
	dates must be non-contiguous, and they are joined by two dots ``'..'``.

	As a last rule, if the start and end date mark the beginning and end of the
	same superordinate period, and if the dates appear to be contiguous, then only this 
	period is returned:
	
	>>> from datetime import datetime as dt, timedelta as td
	>>> dts2str([dt(1986,2,1,0)+i*td(0.25) for i in range(112)])
	'198602'

	Here, the dates span the entire February 1986, but do not extend into March. This 
	example, and all the following assume a time step of 6 hours.

	>>> dts2str([dt(1986,2,1,0)+i*td(0.25) for i in range(236)])
	'198602-198603'

	Here, the dates span the entire February and March of 1986. If only the 
	start and end date are given

	>>> dts2str([dt(1986,2,1,0), dt(1986,3,31,18)])
	'198602..198603'

	the two dates are joined by two dots. The same mechanisms apply also to years

	>>> dts2str([dt(1979,1,1,0)+i*td(0.25) for i in range(51136)])
	'1979-2013'

	and days.

	>>> dts2str([dt(1986,4,7,0)+i*td(0.25) for i in range(4)]
	'19860407'

	However, note that, 

	>>> dts2str([dt(1979,1,1,0)+i*td(0.25) for i in range(51137)])
	'1979-2014010100'

	because 1 January 2014 is the beginning of a new superordinate period rather than the 
	ending of an old. Furthermore,

	>>> dts2str([dt(1986,4,7,6), dt(1986,4,7,18), dt(1986,4,7,0)])
	'19860407..19860407'

	because even if the first and last time step of the superordinate period (which is 7 
	April 1986) is present, not all time steps in between are listed.

	Parameters
	----------
	dates : datetime or set/list/tuple of datetime
	    Dates to be represented by a string. The order of the dates is irrelevant.
	agg : str
	    String representation for the time aggregation interval.
	
	Returns
	-------
	str 
	    Representation of the dates, to be used for example in file names.
	'''

	# Special cases: Only one datetime object given
	if type(dates) == dt:
		return dates.strftime('%Y%m%d%H')
	
	elif len(dates) == 1:
		return dates[0].strftime('%Y%m%d%H')
		
	# General case: Several datetime objects given
	dta = min(dates)
	dtz = max(dates)

	if dta.hour == 0:
		if dta.day == 1:
			if dta.month == 1:
				ret = dta.strftime('%Y')
			else:
				ret = dta.strftime('%Y%m')
		else:
			ret = dta.strftime('%Y%m%d')
	else:
		ret = dta.strftime('%Y%m%d%H')
	
	# Try to find an applicable time aggregation interval, either by argument or from the data set
	if conf.timestep:
		if agg:
			agg = getattr(tagg, agg)
		else: 
			agg = getattr(tagg, conf.timestep)

		agg = agg(dta, dtz)
	else:
		agg = None
	
	# If an aggregation interval is available, find out if the given dates are contiguous 
	# and simplify the date representation accordingly
	if agg:
		ret += '-'
		contiguous = True
		aggdates = list(agg)
		for date in aggdates:
			if not date in dates:
				ret += '..'
				contiguous = False
				break

		last = agg.start_next(aggdates[-1]) - conf.timestep
		if dtz.hour == last.hour:
			if dtz.day == last.day:
				if dtz.month == last.month:
					if dtz.year == dta.year and contiguous:
						ret = dtz.strftime('%Y')
					else:
						ret += dtz.strftime('%Y')
				else:
					if dtz.year == dta.year and dtz.month == dta.month and contiguous:
						ret = dtz.strftime('%Y%m')
					else:
						ret += dtz.strftime('%Y%m')
			else:
				if dtz.year == dta.year and dtz.month == dta.month \
						and dtz.day == dta.day and contiguous:
					ret = dtz.strftime('%Y%m%d')
				else:
					ret += dtz.strftime('%Y%m%d')
		else:
			ret += dtz.strftime('%Y%m%d%H')
	else:
		ret += '..' + dtz.strftime('%Y%m%d%H')

	return ret


# the end
