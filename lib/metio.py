#!/usr/bin/env python
# -*- encoding: utf-8

import os
import sys

import numpy as np
import scipy.io.netcdf as nc3
import Scientific.IO.NetCDF as nc
import scipy.io.matlab as mat
import pytz

from settings import conf as c
import utils
from gridlib import grid_by_nc, grid_by_static

from datetime import datetime as dt, timedelta as td

import dynfor
dynlib_version = (''.join(dynfor.consts.version)).strip()


''' Input/Output library for meteorological data

Provides convenient access to met-data stored in numpy, netCDF and Matlab formats,
and provides functions to save meteorological fields and lines in netCDF files.
'''


# #############################################################################
# 1. Open all sorts of files
# 

# TODO: Restrict the expected file types by providing a file name extension
# TODO: Allow usage without directly requesting a variable, making q optional.
def metopen(filename, q, cut=slice(None), verbose=False, no_dtype_conversion=False, no_static=False):
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
		The requested variable within the file.
	cut : slice
		Obsolete and to be removed
	verbose : bool
		*Optional*, default ``False``. Print debug information on which files are 
		being looked for.
	no_dtype_conversion : bool
		*Optional*, default ``False``. By default, ``metopen`` uncompresses data in the 
		file and converts all data to float64. This behaviour can be suppressed by 
		setting this option to ``True``.
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
		Requested data.
	grid.gridlib
		If ``no_static=False`` meta-information about the requested data.
	'''
	if type(cut) == tuple and len(cut) > 1:
		cuts = tuple(list(cut)[1:])
	else: 
		cuts = slice(None)
	
	tried = []
	for path in c.datapath:
		static = None

		if verbose:
			print 'Trying: '+path+'/'+filename+'.*'

		if os.path.exists(path+'/'+filename+'.npy'):
			dat = np.load(path+'/'+filename+'.npy', mmap_mode='r')
			dat = dat[cut]
			print 'Found '+path+'/'+filename+'.npy'
			f = None
		elif os.path.exists(path+'/'+filename+'.npz'):
			f   = np.load(path+'/'+filename+'.npz')
			if q not in f.files:
				tried.append(path+'/'+filename+'.npz')
				continue
			dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.npz'
		elif os.path.exists(path+'/'+filename+'.mat'):
			f   = mat.loadmat(path+'/'+filename+'.mat')
			if q not in f:
				tried.append(path+'/'+filename+'.mat')
				continue
			dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.mat'
		elif os.path.exists(path+'/'+filename+'.nc'):
			#f   = nc.netcdf_file(path+'/'+filename+'.nc', 'r')
			f   = nc.NetCDFFile(path+'/'+filename+'.nc', 'r')
			var = f.variables[q]
			if q not in f.variables:
				tried.append(path+'/'+filename+'.nc')
				continue
			dat = utils.scale(var, cut)
			if not no_static:
				static = grid_by_nc(f, var, cut=cuts)
				# TODO: Where to search for topography in nc files?
				static.oro = np.zeros((static.ny, static.nx))
				static.oro = static.oro[cuts]
			print 'Found '+path+'/'+filename+'.nc'
		else:
			tried.append(path)
			continue
		
		if not no_dtype_conversion:
			dat = dat.astype('f8')

		if not no_static:
			if not static:
				static = get_static(cuts, verbose, no_dtype_conversion)
						
		else:
			return f, dat

		return f, dat, static
	
	raise RuntimeError, '%s.* not found in any data location. \nTried the following (in order):\n\t%s' % (filename, '\n\t'.join(tried))


# Save dat as a netCDF file, using the metadata in static. 
def metsave(dat, static, time, plev, q, q_units=None, compress_to_short=True):
	''' Save data in a netCDF file

	Parameters
	----------
	dat : np.ndarray
		Data to be saved.
	static : gridlib.grid
		Some meta information about the data, like the grid information.
	time : str
		String representation of the time period covered by the data, e.g. ``'2005'`` for the year
		2005 or ``'20050317'`` for 17 March 2005.
	plev : str
		String representation of the vertical level on which the data is defined, e.g. ``'700'`` 
		for 700 hPa or ``'pv2000'`` for the PV2-surface.
	q : str
		The variable name identifier, following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``.
	q_units : str
		Obsolete.
	compress_to_short : bool
		*Optional*, default ``True``. By default, ``metsave`` compresses the data by converting
		the data field into int16, using the float64 ``add_offset`` and ``scale_factor`` attributes 
		to represent the data.
	'''
	s = dat.shape
	if not len(s) == 3 or not s[1:] == (361,720):
		raise NotImplementedError, 'dat does not seem to be a ERA-Interim like (t,y,x)-array.'
	
	now = dt.now(pytz.timezone('Europe/Oslo'))
	of = nc3.netcdf_file(c.opath+'/'+(c.file_std % {'time': time, 'plev': plev, 'q': c.qi[q]})+'.nc', 'w')
	of._attributes = {'Conventions': 'CF-1.0', 
			'history': '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	}

	of.createDimension('time', dat.shape[0])
	of.createDimension('latitude', dat.shape[1])
	of.createDimension('longitude', dat.shape[2])

	ot = of.createVariable('time', 'i', ('time',))
	ot._attributes = {'long_name': 'time', 'units': static.t_unit}
	ot[:] = static.t
	# TODO: Prescribing the name but taking the units from static is inconsistent!
	olat = of.createVariable('latitude', 'f', ('latitude',))
	olat._attributes = {'long_name': 'latitude', 'units': static.y_unit}
	olat[:] = static.y[:,0]
	olon = of.createVariable('longitude', 'f', ('longitude',))
	olon._attributes = {'long_name': 'longitude', 'units': static.x_unit}
	olon[:] = static.x[0,:]
	
	if compress_to_short:
		ovar = of.createVariable(q, 'h', ('time', 'latitude', 'longitude',))
		dat, scale, off = utils.unscale(dat)
		ovar._attributes = {'long_name': c.q_long[q], 'units': c.q_units[q],
				'add_offset': off, 'scale_factor': scale}
	else:
		ovar = of.createVariable(q, 'f', ('time', 'latitude', 'longitude',))
		ovar._attributes = {'long_name': c.q_long[q], 'units': c.q_units[q]}
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
	of = nc.netcdf_file(c.opath+'/'+(c.file_std % {'time': time, 'plev': plev, 'q': c.qi[q]})+'.nc', 'w')
	of._attributes = {'Conventions': 'CF-1.0', 
			'history': '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	}
	
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
	olon._attributes = {'long_name': static.x_name, 'units': static.x_unit}
	olon[:] = static.x[0,:]
	olat = of.createVariable(static.y_name, 'f', (static.y_name,))
	olat._attributes = {'long_name': static.y_name, 'units': static.y_unit}
	olat[:] = static.y[:,0]

	ot = of.createVariable(static.t_name, 'i', ('time',))
	ot._attributes = {'long_name': static.t_name, 'units': static.t_unit}
	ot[:] = static.t
	olidx = of.createVariable('pointindex', 'i', ('pointindex',))
	olidx._attributes = {'long_name': 'Index of point along all lines', 'units': '1'}
	olidx[:] = range(llen)
	olity = of.createVariable('infotype', 'c', ('infotype',))
	olity._attributes = {'long_name': 'Type of info stored for point', 'units': 'enum'}
	olity[:] = ['X', 'Y', 'I']
	ooidx = of.createVariable('lineindex', 'i', ('lineindex',))
	ooidx._attributes = {'long_name': 'Index of line', 'units': '1'}
	ooidx[:] = range(olen)
	
	oq = of.createVariable(q, 'f', ('time', 'pointindex', 'infotype',))
	oq._attributes = {'long_name': c.q_long[q], 'units': 'mixed'}
	oq[::] = dat[:,:llen,:]

	oqoff = of.createVariable(qoff, 'i', ('time', 'lineindex',))
	oqoff._attributes = {'long_name': 'Index of first point of line', 'units': '1'}
	oqoff[::] = datoff[:,:olen]

	of.close()

	return


# Get static information
def get_static(cuts=slice(None), verbose=False, no_dtype_conversion=False):
	''' Get standard meta-information (for ERA-Interim)

	Parameters
	----------
	cut : slice
		Obsolete and to be removed
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
	fo, oro = metopen('static', 'oro', cuts, verbose, no_dtype_conversion, True)
	static = grid_by_static(fo, cut=cuts)
	static.oro = oro[::]
	fo.close()

	return static


# #############################################################################
# 2. Generalised data fetchers
# 

# TODO: Generalize this function to work on any (constant) input data interval. 
# Or possibly any input data, even with irregular time axis?
def get_instantaneous(q, dates, plevs=None, yidx=None, xidx=None, tavg=False, quiet=False, force=False, **kwargs):
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
	yidx : np.ndarray or slice
		*Optional*. Restrict the output to the given set or range of y-indexes.
	xidx : np.ndarray or slice
		*Optional*. Restrict the output to the given set or range of x-indexes.
	tavg : bool
		*Optional*, default ``False``. Instead of returning a time-dimension, return
		the temporal average of the requested data.
	quiet : bool
		*Optional*, default ``False``. Suppress any output from this function.
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
	# None means "take everything there is as pressure levels"
	if not plevs:
		plevs = c.plevs
	elif not type(plevs) == np.ndarray and not type(plevs) == list:
		plevs = [plevs,]
	
	if type(yidx) == np.ndarray:
		yidxs = yidx
	elif yidx == None:
		yidxs = slice(None)
	else:
		yidxs = yidx
	
	if type(xidx) == np.ndarray:
		xidxs = xidx
	elif xidx == None:
		xidxs = slice(None)
	else:
		xidxs = xidx
	
	dt0 = dt(1979,1,1,0)
	# Convert dates to time indexes
	if type(dates) not in ([np.ndarray, list, tuple, set]):
		dates = [dates, ]
	years = set(map(lambda date: date.year, dates))
	years = np.arange(min(years), max(years)+1)
	tidxs = map(lambda date: date-dt0, dates)
	tidxs = map(lambda diff: diff.days*4 + diff.seconds/21600, tidxs)
	tsmin, tsmax = min(tidxs), max(tidxs)
	
	# Checking max length
	maxtlen = 3000
	if tsmax-tsmin > maxtlen:
		if force:
			print 'Warning: you requested %d time steps.' % (tsmax-tsmin)
		else:
			raise RuntimeError, 'Cowardly refusing to fetch %d > %d time steps.\nUse force=True to override.' % (tsmax-tsmin, maxtlen)
	
	dat = None
	for year in years:
		# Construct the slice
		fst = dt(year,1,1,0) - dt0
		lst = dt(year,12,31,18) - dt0
		fst = fst.days*4 + fst.seconds/21600
		lst = lst.days*4 + lst.seconds/21600 +1
		# Leave out unnecessary indexes for better compatibility
		if not type(xidxs) == np.ndarray and xidxs == slice(None): 
			if yidxs == slice(None):
				cut = (slice(max(tsmin - fst, 0),min(1+tsmax - fst, lst - fst)), )
			else:
				cut = (slice(max(tsmin - fst, 0),min(1+tsmax - fst, lst - fst)), yidxs)
		else:
			cut = (slice(max(tsmin - fst, 0),min(1+tsmax - fst, lst - fst)), yidxs, xidxs)
		datcut = slice(fst+cut[0].start-tsmin, fst+cut[0].stop-tsmin)
		
		# One or more vertical levels?
		i = 0
		for plev in plevs:
			if not quiet:
				print "Reading from "+c.file_std % {'time': year, 'plev': plev, 'q': c.qi[q]}
			if type(dat) == type(None):
				if kwargs.get('no_static', False):
					f, d = metopen(c.file_std % {'time': year, 'plev': plev, 'q': c.qi[q]}, q, cut=cut, **kwargs)
					static = None
				else:
					f, d, static = metopen(c.file_std % {'time': year, 'plev': plev, 'q': c.qi[q]}, q, cut=cut, **kwargs)
					kwargs['no_static'] = True
				s = tuple([1+tsmax-tsmin,len(plevs)] + list(d.shape)[1:])
				dat = np.empty(s, dtype=d.dtype)
			else:
				f, d = metopen(c.file_std % {'time': year, 'plev': plev, 'q': c.qi[q]}, q, cut=cut, **kwargs)
			dat[datcut,i,::] = d
			i += 1

	# Time-averaging if specified
	if tavg and len(dates) > 1:
		dat = dat.mean(axis=0)
	
	dat = dat.squeeze()
	
	return dat, static


# Get aggregated (average, standard deviation, etc.) fields
def get_aggregate(q, year=None, plev=None, yidx=None, xidx=None):
	''' Request aggregate data 
	
	This function is not implemented yet. It is only included to define the future API.
	'''
	raise NotImplementedError, 'If you need it, implement it!'

	return dat


# the end
