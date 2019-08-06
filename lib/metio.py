#!/usr/bin/env python
# -*- encoding: utf-8

''' Input/Output library for meteorological data

Provides convenient access to met-data stored in numpy, netCDF and Matlab formats,
and provides functions to save meteorological fields and lines in netCDF files.
'''

from __future__ import absolute_import, unicode_literals, print_function

import os
import sys
import glob

import numpy as np
import netCDF4 as nc
import scipy.io.matlab as mat
import pytz
import calendar
from datetime import datetime as dt, timedelta as td
from dateutil.relativedelta import relativedelta as rtd

from . import settings_basic as se
from . import utils
from .gridlib import grid_by_nc, grid_by_static
from . import tagg
from . import dynfor

from . import version as dynlib_version



# Maxmimum number of time steps to be requested at once in get_instantaneous 
# and get_aggregate without having to use force
MAX_TLEN = 3000
NO_ENDING = 'nc' # ETH compatibility: if the file exists without file ending, try to interpret it as a netCDF file


# #############################################################################
# 1. Open all sorts of files
# 

# TODO: Restrict the expected file types by providing a file name extension
def metopen(filename, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False, no_static=False, quiet=False, mode='r'):
	''' Find and open files by name
	
	Uses the se.conf.datapath list to locate files in a variety of different locations. 
	In user scripts and user settings files, this variable will typically be available
	via conf.datapath.

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
	mode :
	    *Optional*, default ``'r'``. Only effective for netCDF files. The read/write mode
	    with which to open the file. Valid values are ``'r'`` for read-only access, `'a'``
	    or ``'r+'` for read-write access and and ``'w'``for replacing the given file.
	
	Returns
	-------
	data file object
	    python data, netCDF or Matlab file object.
	np.ndarray
	    If q given, data of the requested variable.
	grid.gridlib
	    If ``no_static=False`` meta-information about the requested data.
	'''

	def handle_npy(filepath):
		if not mode == 'r':
			print('WARNING: Can only open npy files in read mode!')
		if q:
			dat = np.load(filepath, mmap_mode='r', encoding='bytes')
			dat = dat[cut]
		else:
			dat = None
		if not quiet:
			print('Found '+filepath)
		f = None

		return f, dat

	def handle_npz(filepath):
		if not mode == 'r' and not quiet:
			print('WARNING: Can only open npz files in read mode!')
		f = np.load(filepath, encoding='bytes')
		if q:
			if q not in f.files:
				tried.append(filepath)
			dat = f[q][cut]
		else:
			dat = None
		if not quiet:
			print('Found '+filepath)

		return f, dat
	
	def handle_mat(filepath):
		if not mode == 'r' and not quiet:
			print('WARNING: Can only open mat files in read mode!')
		f = mat.loadmat(filepath)
		if q:
			if q not in f:
				tried.append(filepath)
			dat = f[q][cut]
		else:
			dat = None
		if not quiet:
			print('Found '+filepath)

		return f, dat

	def handle_nc(filepath):
		f = nc.Dataset(filepath, mode)
		f.set_auto_scale(False)
		f.set_auto_mask(False)
		if q:
			var = f.variables[q]
			if q not in f.variables:
				tried.append(filepath)
			if not no_dtype_conversion:
				dat = utils.scale(var, cut=cut)
			else:
				dat = var[cut]
		else:
			dat = None

		if not no_static:
			if q:
				static = grid_by_nc(f, var)
			else:
				static = grid_by_nc(f)
			# TODO: Where to search for topography in nc files?
			static.oro = np.zeros((static.ny, static.nx))
		else:
			static = None

		if not quiet:
			print('Found '+filepath)

		return f, dat, static

	if not type(cut) == slice:
		raise ValueError('cut must be a 1-dimensional slice object')
	
	tried = []
	for path in se.conf.datapath:
		static = None

		if verbose:
			print('Trying: '+path+'/'+filename+'.*')

		if os.path.exists(path+'/'+filename+'.npy'):
			f, dat = handle_npy(path+'/'+filename+'.npy')
		elif os.path.exists(path+'/'+filename) and NO_ENDING == 'npy':
			f, dat = handle_npy(path+'/'+filename)

		elif os.path.exists(path+'/'+filename+'.npz'):
			f, dat = handle_npz(path+'/'+filename+'.npz')
		elif os.path.exists(path+'/'+filename) and NO_ENDING == 'npz':
			f, dat = handle_npz(path+'/'+filename)

		elif os.path.exists(path+'/'+filename+'.mat'):
			f, dat = handle_mat(path+'/'+filename+'.mat')
		elif os.path.exists(path+'/'+filename) and NO_ENDING == 'mat':
			f, dat = handle_mat(path+'/'+filename)

		elif os.path.exists(path+'/'+filename+'.nc'):
			f, dat, static = handle_nc(path+'/'+filename+'.nc')
		elif os.path.exists(path+'/'+filename) and NO_ENDING == 'nc':
			f, dat, static = handle_nc(path+'/'+filename)

		else:
			tried.append(path)
			continue
		
		if q and not no_dtype_conversion and not dat.dtype == 'f8':
			dat = dat.astype('f8')

		if not no_static:
			if not static:
				static = get_static(verbose, no_dtype_conversion, quiet)
		else:
			if q:
				return f, dat
			else:
				return f
		if q:
			return f, dat, static
		else:
			return f, static
	
	raise ValueError('%s.* not found in any data location. \n'
			'Tried the following (in order):\n\t%s' % (filename, '\n\t'.join(tried)) )


def metdiscover(filename):
	''' Find files by name, potentially including wildcards or other patterns
	
	Parameters
	----------
	filename : str
	    The file name (pattern) to be found.

	Returns
	-------
	List of str
	    A list of filenames that match the given filename pattern. Might be empty.
	'''

	found = set([])
	for path in se.conf.datapath:
		_found = glob.glob(path+'/'+filename+'.nc')
		found.update([os.path.basename(_f[:-3]) for _f in _found])

		_found = glob.glob(path+'/'+filename+'.npy')
		_found.extend(glob.glob(path+'/'+filename+'.npz'))
		_found.extend(glob.glob(path+'/'+filename+'.mat'))
		found.update([os.path.basename(_f[:-4]) for _f in _found])

	return found


# Save dat as a netCDF file, using the metadata in static. 
def metsave(dat, static, q, plev, agg=None, compress_to_short=True):
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
	agg : str
	    *Optional*. String representation for the time aggregation interval, e.g. ``'daily'`` or 
	    ``'cal_month'``. Affects the file name, and the string representation of the given dates.
	compress_to_short : bool
	    *Optional*, default ``True``. By default, ``metsave`` compresses the data by converting
	    the data field into int16, using the float64 ``add_offset`` and ``scale_factor`` attributes 
	    to represent the data.
	'''
	s = dat.shape

	# TODO: Why check against gridsize? Or at least: Why not allowing gridsize to be unset?
	if not plev == 'sfc' and (not len(s) == 4 or (se.conf.gridsize and not s[2:] == se.conf.gridsize)):
		raise NotImplementedError('dat does not seem to be a context-conform (t,z,y,x)-array.'
				'\n Expected 4D with horizontal grid of dimension %s'
				'\n Got shape: %s' % (se.conf.gridsize, s) )
	elif plev == 'sfc' and (not len(s) == 3 or (se.conf.gridsize and not s[1:] == se.conf.gridsize)):
		raise NotImplementedError('dat does not seem to be a context-conform surface (t,y,x)-array.'
				'\n Expected 3D with horizontal grid of dimension %s'
				'\n Got shape: %s' % (se.conf.gridsize, s) )
	
	if not se.conf.oformat == 'nc':
		raise NotImplementedError('Currently only saving in netCDF implemented in metsave.')

	if not s[0] == len(static.t):
		raise ValueError('Time dimension in data (%s) and static (%s) are not equally long.' % (s[0], len(static.t)))
	if not plev == 'sfc' and not s[1] == len(static.z):
		raise ValueError('z-dimension in data (%s) and static (%s) are not equally long.' % (s[1], len(static.z)))
	if not s[-2] == len(static.y[:,0]):
		raise ValueError('y-dimension in data (%s) and static (%s) are not equally long.' % (s[-2], len(static.y[:,0])))
	if not s[-1] == len(static.x[0,:]):
		raise ValueError('x-dimension in data (%s) and static (%s) are not equally long.' % (s[-1], len(static.x[0,:])))

	now = dt.now(pytz.timezone(se.conf.local_timezone))
	if not agg:
		of = nc.Dataset(se.conf.opath+'/'+(se.conf.file_std % {
				'time': dts2str(static.t_parsed), 'plev': plev, 'qf': se.conf.qf[q]})+'.nc', 
			'w', format='NETCDF3_CLASSIC')
	else:
		of = nc.Dataset(se.conf.opath+'/'+(se.conf.file_agg % { 'agg': agg,
				'time': dts2str(static.t_parsed, agg), 'plev': plev, 'qf': se.conf.qf[q]})+'.nc', 
			'w', format='NETCDF3_CLASSIC')
	of.setncatts({'Conventions': 'CF-1.0', 
			'history': '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	})

	of.createDimension('time', None)
	if not plev == 'sfc':
		of.createDimension(static.z_name, s[1])
	
	if getattr(static, 'rotated', False):
		y_name = static.rot_y_name
		x_name = static.rot_x_name
	else:
		y_name = static.y_name
		x_name = static.x_name
	
	of.createDimension(y_name, s[-2])
	of.createDimension(x_name, s[-1])

	if not plev == 'sfc':
		known_vertical_level_units = {
				'Pa': ('pressure', 'down'),
				'K': ('isentropic', 'up'),
				'PVU': ('potential_vorticity', 'up'),
		}
		if not static.z_unit in known_vertical_level_units:
			raise ValueError('Unknown vertical level type unit: `%s`' % static.z_unit)
		z_name, z_positive = known_vertical_level_units[static.z_unit]

	ot = of.createVariable('time', 'i', ('time',))
	ot.setncatts({'long_name': 'time', 'units': static.t_unit})
	ot[::] = static.t
	if not plev == 'sfc':
		olev = of.createVariable(static.z_name, 'f', (static.z_name,))
		olev.setncatts({'long_name': z_name, 'units': static.z_unit, 'axis': 'Z', 'positive': z_positive})
		olev[::] = static.z[:]
		dims = ('time', static.z_name, y_name, x_name,)
	else:
		dims = ('time', y_name, x_name,)
	
	if getattr(static, 'rotated', False):
		olat = of.createVariable(static.rot_y_name, 'f', (y_name,))
		olat.setncatts({'long_name': static.rot_y_longname, 'units': static.y_unit, 'axis': 'Y'})
		olat[::] = static.rot_y[:,0]
		olon = of.createVariable(static.rot_x_name, 'f', (x_name,))
		olon.setncatts({'long_name': static.rot_x_longname, 'units': static.x_unit, 'axis': 'X'})
		olon[::] = static.rot_x[0,:]

		olat = of.createVariable(static.y_name, 'f', (y_name, x_name,))
		olat.setncatts({'long_name': static.y_name, 'units': static.y_unit})
		olat[::] = static.y
		olon = of.createVariable(static.x_name, 'f', (y_name, x_name,))
		olon.setncatts({'long_name': static.x_name, 'units': static.x_unit})
		olon[::] = static.x

		orot = of.createVariable('rotated_pole', 'i', ())
		orot.setncatts({
			'grid_north_pole_longitude': static.rot_np[1], 
			'grid_north_pole_latitude': static.rot_np[0],
		})

	else:
		olat = of.createVariable(static.y_name, 'f', (y_name,))
		olat.setncatts({'long_name': static.y_name, 'units': static.y_unit, 'axis': 'Y'})
		olat[::] = static.y[:,0]
		olon = of.createVariable(static.x_name, 'f', (x_name,))
		olon.setncatts({'long_name': static.x_name, 'units': static.x_unit, 'axis': 'X'})
		olon[::] = static.x[0,:]
	
	if compress_to_short:
		dat, scale, off, fill = utils.unscale(dat)
		if fill: 
			ovar = of.createVariable(q, 'i2', dims, fill_value=fill)
			ovar.set_auto_scale(False)

			ovar.missing_value = fill
		else:
			ovar = of.createVariable(q, 'i2', dims)

		ovar.setncatts({'add_offset': off, 'scale_factor': scale})
	else:
		ovar = of.createVariable(q, dat.dtype, dims)
	
	ovar.setncatts({'long_name': se.conf.q_long[q], 'units': se.conf.q_units[q]})
	# Add some attributes to make ncview display the rotated grid correctly
	if getattr(static, 'rotated', False):
		ovar.setncatts({'grid_mapping': 'rotated_pole', 'coordinates': static.x_name+' '+static.y_name})
	ovar[::] = dat

	of.close()

	return


# Save dat as a netCDF file, using the metadata in static. 
def metsave_lines(dat, datoff, static, plev, q, qoff, additional_axes=[]):
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
	plev : str
	    String representation of the vertical level on which the data is defined, e.g. ``'700'`` 
	    for 700 hPa or ``'pv2000'`` for the PV2-surface.
	q : str
	    The variable name identifier, e.g. ``'jetaxis'``.
	qoff : str
	    The variable name identifier for the list of offsets, e.g. ``'jaoff'``. 
	additional_axes : list of dict
	    List of additional axes, in order, to be placed between time and space.
	    The dict must have the keys ``'name'`` and ``'data'``, defining the axis name and its
	    data values. All additional keys will be taken over as attributes for the axis variable.
	'''
	if not len(dat.shape) == (3+len(additional_axes)) and len(datoff.shape) == (2+len(additional_axes)):
		raise RuntimeError('dat and/or datoff have the wrong number of dimensions')
	if not dat.shape[:-2] == datoff.shape[:-1]:
		raise RuntimeError('dat and datoff have different axis lengths')
	if not dat.shape[-1] == 3:
		raise RuntimeError('dat does not have size 3 in the last dimension')

	now = dt.now(pytz.timezone('Europe/Oslo'))
	filename = se.conf.opath+'/'+(se.conf.file_std % {'time': dts2str(static.t_parsed), 
		'plev': plev, 'qf': se.conf.qf[q]})+'.nc'
	print('Saving to %s' % filename)
	of = nc.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
	of.setncatts({'Conventions': 'CF-1.0', 
			'history': '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	})
	
	llen = dat.shape[-2]
	olen = datoff.shape[-1]

	of.createDimension(static.t_name, dat.shape[0])
	of.createDimension('pointindex', llen)
	of.createDimension('infotype', 3)
	of.createDimension('lineindex', olen)
	
	# Define helper dimensions; not explicitly used in any variable but useful to interpret the grid point indexes stored in dat
	of.createDimension(static.x_name, static.x.shape[1])
	of.createDimension(static.y_name, static.y.shape[0])
	olat = of.createVariable(static.y_name, 'f', (static.y_name,))
	olat.setncatts({'long_name': static.y_name, 'units': static.y_unit})
	olat[:] = static.y[:,0]
	olon = of.createVariable(static.x_name, 'f', (static.x_name,))
	olon.setncatts({'long_name': static.x_name, 'units': static.x_unit})
	olon[:] = static.x[0,:]
	
	add_dimnames = ()
	for dim in additional_axes:
		if 'name' not in dim or 'data' not in dim:
			raise ValueError('Additional dimensions must have a name and data')
		of.createDimension(dim['name'], len(dim['data']))
		odim = of.createVariable(dim['name'], 'f', (dim['name'], ))
		odim[:] = dim['data']
		add_dimnames += (dim['name'], ) 

		del dim['name']
		del dim['data']
		odim.setncatts(dim)

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
	
	oq = of.createVariable(q, 'f', ('time', ) + add_dimnames + ('pointindex', 'infotype',))
	oq.setncatts({'long_name': se.conf.q_long[q], 'units': 'mixed'})
	oq[::] = dat[::]

	oqoff = of.createVariable(qoff, 'i', ('time', ) + add_dimnames + ('lineindex',))
	oqoff.setncatts({'long_name': 'Index of first point of line', 'units': '1'})
	oqoff[::] = datoff[::]

	of.close()

	return


def metsave_timeless(dat, static, name, ids=None, q=None, plev=None, compress_to_short=True, 
		timeseries={}, global_atts={}):
	''' Save time-independent data like composites or EOFs in a netCDF file

	The data is saved either to an existing file with a matching name in se.conf.datapath, 
	or, if such a file does not exist, to a new file in se.conf.opath.

	Optionally, time information may be retained in time series. Time series may use the ID
	dimension, if ids set.

	Parameters
	----------
	dat : (dict of) np.ndarray 
	    Data to be saved. If an numpy array, the variable name q used for meta information.
	static : gridlib.grid
	    Some meta information about the data, like the grid information.
	name : str
	    Name of the collection, identifying for example the composite or EOF.
	ids : list of str
	    *Optional*. Identifying names for each index in the ``id`` dimension.
	q : str
	    *Optional*. Only used and required if dat is a numpy array. The variable name identifier, 
	    following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``.
	plev : str
	    *Optional*. Only used and required if dat is a numpy array. String representation of the 
	    vertical level on which the data is defined, e.g. ``'700'`` for 700 hPa or ``'pv2000'`` 
	    for the PV2-surface.
	compress_to_short : bool
	    *Optional*, default ``True``. By default, ``metsave`` compresses the data by converting
	    the data field into int16, using the float64 ``add_offset`` and ``scale_factor`` attributes 
	    to represent the data.
	'''

	if not type(dat) == dict:
		if not q or not plev:
			raise ValueError('Variable name and vertical level required, if dat is not a dict!')
		datdict = {(plev, q): dat}
	else:
		datdict = dat
	
	s = static.x.shape
	if ids:
		s = (len(ids), ) + s
	
	now = dt.now(pytz.timezone(se.conf.local_timezone))
	history = '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	filename = se.conf.file_timeless % {'time': dts2str(static.t_parsed), 'name': name}
	try:
		f = metopen(filename, no_static=True, mode='a')
		print('Saving to existing %s' % filename)
		new = False

	except ValueError: 
		f = nc.Dataset(se.conf.opath+'/'+filename+'.nc', 'w', format='NETCDF4')
		print('Saving to %s/%s.nc' % (se.conf.opath, filename))
		new = True
	
	# Consistency checks
	if not new:
		for att in global_atts:
			exist = getattr(f, att, None)
			if exist:
				if not exist == global_atts[att]:
					raise ValueError('Existing file found with incompatible attributes: '
							'`%s` (existing: %s, given: %s)' % (
								att, exist, global_atts[att]))
			else:
				setattr(f, att, global_atts[att])

		if timeseries and not 'time' in f.dimensions:
			raise ValueError('Time series given, but no time dimension present in existing file.')
		
		if ids and 'id' in f.variables:
			for id1, id2 in zip(f.variables['id'], ids):
				if not id1 == id2:
					raise ValueError('Existing file found with incompatible ID dimension (existing: %s, given: %s)' % (f.variables['id'][::], ids))
			data_dimensions = ('id', )
		elif ids and not 'id' in f.variables:
			raise ValueError('Existing file does not contain an ID dimension')
		elif not ids and 'id' in f.variables:
			raise ValueError('Existing file contains ID dimension, but no ids given!')
		else:
			data_dimensions = ()

		if not static.y_name in f.variables:
			raise ValueError('Dimension `%s` not found in existing file.' % static.y_name)
		if not static.x_name in f.variables:
			raise ValueError('Dimension `%s` not found in existing file.' % static.x_name)

		data_dimensions += (static.y_name, static.x_name, )
	
	# Create netCDF dimensions and attributes
	else:
		f.setncatts(global_atts)
		f.setncatts({'Conventions': 'CF-1.0', 
				'history': history, 
		})

		if timeseries:
			f.createDimension('time', len(static.t))
			ot = f.createVariable('time', 'i', ('time',))
			ot.setncatts({'long_name': 'time', 'units': static.t_unit})
			ot[::] = static.t
		
		if ids:
			f.createDimension('id', s[0])
			f.createDimension(static.y_name, s[1])
			f.createDimension(static.x_name, s[2])

			oid = f.createVariable('id', str, ('id',))
			oid.setncatts({'long_name': 'Identifying name'})
			oid[::] = np.array(ids)

			data_dimensions = ('id', static.y_name, static.x_name, )
		else:
			f.createDimension(static.x_name, s[1])
			f.createDimension(static.y_name, s[0])

			data_dimensions = (static.y_name, static.x_name, )

		olat = f.createVariable(static.y_name, 'f', (static.y_name,))
		olat.setncatts({'long_name': static.y_name, 'units': static.y_unit})
		olat[::] = static.y[:,0]
		olon = f.createVariable(static.x_name, 'f', (static.x_name,))
		olon.setncatts({'long_name': static.x_name, 'units': static.x_unit})
		olon[::] = static.x[0,:]

	for plev, q_ in datdict:
		# Finding 
		if plev:
			if plev in f.groups and q_ in f.groups[plev].variables:
				print('Warning: variable /%s/%s already present, skipping!' % (plev, q_))
				continue
	
			ncvarname = '/%s/%s' % (plev, q_)

		else:
			if q_ in f.variables:
				print('Warning: variable %s already present, skipping!' % q_)
				continue

			ncvarname = q_
		
		compress_to_short_ = compress_to_short
		if q_[-5:] == '_mean':
			q = q_[:-5]
			prefix = 'Mean '
		elif q_[-4:] == '_cnt':
			q = q_[:-4]
			prefix = 'Number in composite of '
			compress_to_short_ = False
		elif q_[-4:] == '_std':
			q = q_[:-4]
			prefix = 'Standard deviation of '
		elif q_[-5:] == '_hist':
			q = q_[:-5]
			prefix = 'Histogram of '
			compress_to_short_ = False
		elif q_[-4:] == '_mfv':
			q = q_[:-4]
			prefix = 'Most frequent value of '
		elif q_[-4:] == '_min':
			q = q_[:-4]
			prefix = 'Minimum '
		elif q_[-4:] == '_max':
			q = q_[:-4]
			prefix = 'Maximum '
		elif q_[-8:] == '_pattern':
			q = q_[:-8]
			prefix = 'EOF loading pattern for '
		else:
			q = q_
			prefix = ''
		
		s_ = datdict[plev,q_].shape
		if not len(s_) == len(data_dimensions):
			raise ValueError('Data for (%s,%s) does not have the required number of dimensions. '
					'Expected: %d, got: %d.' % (
					plev, q_, len(data_dimensions), len(s_)) )
		data_dimensions_ = ()
		squeezeme = False
		for dimname, dim in zip(data_dimensions, s_):
			if dim > 1:
				data_dimensions_ += (dimname,)
			else:
				squeezeme = True

		if squeezeme:
			dat_ = datdict[plev,q_].squeeze()
			compress_to_short_ = False
		else: 
			dat_ = datdict[plev,q_]

		if compress_to_short_: 
			dat_, scale, off, fill = utils.unscale(dat_)
			if fill: 
				ovar = f.createVariable(ncvarname, 'i2', data_dimensions_, fill_value=fill)
				ovar.set_auto_scale(False)

				ovar.missing_value = fill
			else:
				ovar = f.createVariable(ncvarname, 'i2', data_dimensions_)

			ovar.setncatts({'long_name': prefix+se.conf.q_long[q], 'units': se.conf.q_units[q], 'history': history,
					'add_offset': off, 'scale_factor': scale})
		else:
			ovar = f.createVariable(ncvarname, dat_.dtype, data_dimensions_)
			ovar.setncatts({'long_name': prefix+se.conf.q_long[q], 'units': se.conf.q_units[q], 'history': history})
	
		ovar[::] = dat_

	for plev, q in timeseries:
		# Finding 
		if plev:
			if plev in f.groups and q in f.groups[plev].variables:
				print('Warning: variable /%s/%s already present, skipping!' % (plev, q))
				continue
	
			ncvarname = '/%s/%s' % (plev, q)

		else:
			if q in f.variables:
				print('Warning: variable %s already present, skipping!' % q)
				continue

			ncvarname = q

		if len(timeseries[plev,q].shape) == 2:
			data_dimensions = ('time', 'id')
		else:
			data_dimensions = ('time', )
		
		ovar = f.createVariable(ncvarname, 'f', data_dimensions)
		ovar.setncatts({'long_name': se.conf.q_long[q], 'units': se.conf.q_units[q], 'history': history})
		ovar[::] = timeseries[plev,q]

	f.close()

	return



# Get static information
def get_static(verbose=False, no_dtype_conversion=False, quiet=False):
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
	gridlib.grid
	    Some meta information about the data, like the grid information.
	'''
	
	if not se.conf.file_static:
		raise ValueError('Static file must be configured for get_static() to work!')

	fo, oro = metopen(se.conf.file_static, 'oro', verbose=verbose, no_dtype_conversion=no_dtype_conversion, quiet=quiet, no_static=True)
	static = grid_by_static(fo)
	static.oro = oro[::]
	fo.close()

	return static


# #############################################################################
# 2. Generalised data fetchers
# 

def get_instantaneous(q, dates, plevs=None, tavg=False, force=False, **kwargs):
	''' Data fetcher for instantaneous or periodically averaged fields

	Allows general data requests in the configured data base, e.g. ERA-Interim. The request
	can span several files, e.g. by including several vertical levels or by covering several
	years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
	
	The method internally uses metopen to locate data files. Hence, it will find data in the 
	locations given se.conf.datapath (in user scripts and user settings files typically 
	available as conf.datapath).

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
	    for the PV2-surface. The parameter is only effective if data is split into seperate 
	    files per vertical level.
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
	
	if type(plevs) == type(None):
		plevs = se.conf.plevs
	elif not type(plevs) == np.ndarray and not type(plevs) == list:
		plevs = [plevs,]
	
	if se.conf.years:
		dt0 = dt(se.conf.years[0],1,1,0)
	elif se.conf.times: 
		dt0 = se.conf.times[0]
	else:
		raise ValueError('Could not determine start date of the data set, either se.conf.years or se.conf.times must be provided.')

	dts = se.conf.timestep.total_seconds()

	# Convert dates to time indexes
	if type(dates) not in ([np.ndarray, list, tuple, set]):
		dates = [dates, ]
	years = set(map(lambda date: date.year, dates))
	years = np.arange(min(years), max(years)+1)
	tidxs = [int((date-dt0).total_seconds()/dts) for date in dates]
	tsmin, tsmax = min(tidxs), max(tidxs)

	# Checking max length
	if tsmax-tsmin > MAX_TLEN:
		if force:
			print('Warning: you requested %d time steps.' % (tsmax-tsmin))
		else:
			raise ValueError('Cowardly refusing to fetch %d > %d time steps.\nUse force=True to override.' % (tsmax-tsmin, MAX_TLEN))
	
	# Remove no_static if present
	kwargs.pop('no_static', None)
	static = None

	dat = None
	for year in years:
		# Construct the slice
		fst = dt(year,1,1,0) - dt0
		lst = (dt(year+1,1,1,0) - se.conf.timestep) - dt0
		fst = int(fst.total_seconds()/dts)
		lst = int(lst.total_seconds()/dts) + 1

		# Leave out unnecessary indexes for better compatibility
		cut = slice(max(tsmin - fst, 0),min(1+tsmax - fst, lst - fst))
		datcut = slice(fst+cut.start-tsmin, fst+cut.stop-tsmin)
		
		# One or more vertical levels?
		i = 0
		if type(dat) == type(None):
			f, d, static = metopen(se.conf.file_std % {'time': year, 'plev': plevs[0], 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
			# Inject meta data for npy files
			if not f:
				static.t = np.arange(d.shape[0])*dts/3600
				static.t_parsed = [dt(year,1,1)+i*se.conf.timestep for i in range(d.shape[0])]
			if len(d.shape) == 4 and d.shape[1] > 1:
				separate_plevs = False
				s = (1+tsmax-tsmin, ) + d.shape[1:]
			elif len(d.shape) == 4: 
				separate_plevs = True
				s = (1+tsmax-tsmin, len(plevs), ) + d.shape[2:]
			else:
				separate_plevs = True
				s = (1+tsmax-tsmin, len(plevs), ) + d.shape[1:]

			dat = np.empty(s, dtype=d.dtype)
			if separate_plevs:
				if len(d.shape) < 4:
					d = d[:,np.newaxis,::]
				dat[datcut,0:1,::] = d
			else:
				dat[datcut,::] = d

			static.t = static.t[cut]
			if type(static.t_parsed) == np.ndarray:
				static.t_parsed = static.t_parsed[cut]

			i = 1
		
		if separate_plevs:
			for plev in plevs[i:]:
				f, d, static_ = metopen(se.conf.file_std % {'time': year, 'plev': plev, 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
				if len(d.shape) < 4:
					d = d[:,np.newaxis,::]
				dat[datcut,i:i+1,::] = d
				if i == 0:
					static.t = np.concatenate((static.t, static_.t[cut]))
					if type(static.t_parsed) == np.ndarray:
						static.t_parsed = np.concatenate((static.t_parsed, static_.t_parsed[cut]))
				elif year == years[0]:
					if not static.z_name == static_.z_name or \
							not static.z_unit == static_.z_unit:
						print('WARNING: concatenating vertical levels of different types!')
						static.z_unit = u'MIXED!'
					static.z = np.concatenate((static.z, static_.z))
				i += 1

		elif i == 0:
			f, dat[datcut,::], static_ = metopen(se.conf.file_std % {'time': year, 'plev': plevs[0], 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
			static.t = np.concatenate((static.t, static_.t[cut]))
			if type(static.t_parsed) == np.ndarray:
				static.t_parsed = np.concatenate((static.t_parsed, static_.t_parsed[cut]))

	# Time-averaging if specified
	if tavg and len(dates) > 1:
		dat = dat.mean(axis=0)
	
	#dat = dat.squeeze()
	
	return dat, static


# Get aggregated (average, standard deviation, etc.) fields
def get_aggregate(q, dates, agg, plevs=None, tavg=False, force=False, **kwargs):
	''' Data fetcher for perodically time-averaged fields

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
	agg : str
	    String representation for the time aggregation interval, e.g. ``'daily'`` or 
	    ``'cal_month'``.
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
		plevs = se.conf.plevs
	elif not type(plevs) == np.ndarray and not type(plevs) == list:
		plevs = [plevs,]
	
	if se.conf.years:
		dt0 = dt(se.conf.years[0],1,1,0)
	elif se.conf.times: 
		dt0 = se.conf.times[0]
	else:
		raise ValueError('Could not determine start date of the data set, either se.conf.years or se.conf.times must be provided.')
	
	# Convert dates to time indexes
	if type(dates) not in ([np.ndarray, list, tuple, set]):
		dates = [dates, ]
	years = set(map(lambda date: date.year, dates))
	years = np.arange(min(years), max(years)+1)
	t_iter = getattr(tagg, agg)(dt0)
	t_iter.dtend = t_iter.start_next(max(dates))
	t_list = list(t_iter)

	t_fst = t_iter.start(min(dates))
	t_lst = t_iter.end(max(dates))
	tsmin = t_list.index(t_fst)
	tsmax = t_list.index(t_lst) + 1

	# Checking max length
	if tsmax-tsmin > MAX_TLEN:
		if force:
			print('Warning: you requested %d time steps.' % (tsmax-tsmin))
		else:
			raise ValueError('Cowardly refusing to fetch %d > %d time steps.\nUse force=True to override.' % (tsmax-tsmin, MAX_TLEN))
	
	# Find the position of the year information from a test file name
	MARKER = '__DYNLIB_T_MARKER__'
	testname = se.conf.file_agg % {'agg': agg, 'time': MARKER, 'plev': plevs[0], 'qf': se.conf.qf[q]}
	idx = testname.index(MARKER)

	# Find aggregate file names (multi-year files) that are available on all the requested levels
	filenames = metdiscover(se.conf.file_agg % {'agg': agg, 'time': '[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9]', 
			'plev': plevs[0], 'qf': se.conf.qf[q]})
	for plev in plevs[1:]:
		filenames = filenames.intersection(
			metdiscover(se.conf.file_agg % {'agg': agg, 'time': '[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9]', 
					'plev': plev, 'qf': se.conf.qf[q]})
		)
	
	avail_times = [f[idx:idx+9] for f in filenames]
	avail_years = {}
	for avail_time in avail_times:
		for yr in range(int(avail_time[:4]), int(avail_time[5:9])+1):
			if not yr in avail_years:
				avail_years[yr] = avail_time
	
	# Find aggregate file names (single-year files) that are available on all the requested levels
	filenames2 = metdiscover(se.conf.file_agg % {'agg': agg, 'time': '[0-9][0-9][0-9][0-9]', 
			'plev': plevs[0], 'qf': se.conf.qf[q]})
	for plev in plevs[1:]:
		filenames2 = filenames2.intersection(
			metdiscover(se.conf.file_agg % {'agg': agg, 'time': '[0-9][0-9][0-9][0-9]', 
					'plev': plev, 'qf': se.conf.qf[q]})
		)
	
	for yr in [int(f[idx:idx+4]) for f in filenames2]:
		if not yr in avail_years:
			avail_years[yr] = f[idx:idx+4]
	
	# Check if every requested year is covered
	for yr in years:
		if yr not in avail_years:
			raise ValueError('Could not find data for year %d' % yr)
	
	# Seems like we can finally start reading some data
	
	# Remove no_static if present
	kwargs.pop('no_static', None)
	static = None

	dat = None
	for year in years:
		# Construct the slice
		fst = t_iter.start_after_or_on(dt(year,1,1,0))
		fst = t_list.index(fst)
		lst = t_iter.end_after(dt(year+1,1,1,0)) - t_iter.interval
		if lst > t_list[-1]:
			t_iter.cur = t_list[0]
			t_iter.dtend = lst + t_iter.interval
			t_list = list(t_iter) 
		lst = t_list.index(lst) + 1

		# Calculate offset for several years in one file
		time = avail_years[year]
		fst_file = t_iter.start_after_or_on(dt(int(time[:4]),1,1,0))
		fst_file = t_list.index(fst_file)
		offset = fst - fst_file

		# Leave out unnecessary indexes for better compatibility
		cut = slice(max(tsmin - fst, 0)+offset, min(tsmax - fst, lst - fst)+offset)
		datcut = slice(fst+cut.start-tsmin-offset, fst+cut.stop-tsmin-offset)

		# One or more vertical levels?
		i = 0
		for plev in plevs:
			if type(dat) == type(None):
				f, d, static = metopen(se.conf.file_agg % {'agg': agg, 'time': time, 'plev': plev, 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
				s = (tsmax-tsmin, len(plevs), ) + d.shape[1:]
				dat = np.empty(s, dtype=d.dtype)
				static.t = static.t[cut]
				if type(static.t_parsed) == np.ndarray:
					static.t_parsed = static.t_parsed[cut]
			else:
				f, d, static_ = metopen(se.conf.file_agg % {'agg': agg, 'time': time, 'plev': plev, 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
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
	
	The calculation of the last time step depends on the ``se.conf.timestep`` property.
	For example the last time step of October 2014 in the ERA-Interim data is 18 UTC
	on 31 October 2014.

	If the given number of dates given matches the number of dates expected for the given 
	start and end dates with the ``se.conf.timestep``, the dates are checked wether they are
	contiguous. If yes, the two dates are joined by a minus sign ``'-'``. Otherwise, for a
	non-contiguous set of dates, the start and end dates are joined by two dots ``'..'``.

	As a last rule, if the start and end date mark the beginning and end of the
	same superordinate period, and if the dates appear to be contiguous, then only this 
	period is returned:
	
	>>> from datetime import datetime as dt, timedelta as td
	>>> dts2str([dt(1986,2,1,0)+i*td(0.25) for i in range(112)], 'six_hourly')
	'198602'

	Here, the dates span the entire February 1986, but do not extend into March. This 
	example, and all the following assume a time step of 6 hours.

	>>> dts2str([dt(1986,2,1,0)+i*td(0.25) for i in range(236)], 'six_hourly')
	'198602-198603'

	Here, the dates span the entire February and March of 1986. If only the 
	start and end date are given

	>>> dts2str([dt(1986,2,1,0), dt(1986,3,31,18)], 'six_hourly')
	'198602..198603'

	the two dates are joined by two dots. The same mechanisms apply also to years

	>>> dts2str([dt(1979,1,1,0)+i*td(0.25) for i in range(2924)], 'six_hourly')
	'1979-1980'

	and days.

	>>> dts2str([dt(1986,4,7,0)+i*td(0.25) for i in range(4)], 'six_hourly')
	'19860407'

	However, note that, 

	>>> dts2str([dt(1979,1,1,0)+i*td(0.25) for i in range(2925)], 'six_hourly')
	'1979010100-1981010100'

	because 1 January 2014 is the beginning of a new superordinate period rather than the 
	ending of an old. Furthermore,

	>>> dts2str([dt(1986,4,7,6), dt(1986,4,7,18), dt(1986,4,7,0)], 'six_hourly')
	'19860407..19860407'

	because even if the first and last time step of the superordinate period (which is 7 
	April 1986) is present, not all time steps in between are listed.

	The time step can be overwritten by using the optional argument agg to specify a time
	aggregation interval. In most cases it can be omitted. In this case, the aggregation
	period is set to the time step configured in the context.

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

	def date_format(date, repr_interval):
		''' Unfortunately datetime.strftime does not support dates before 1900, so
		need to implement the required formatter myself. 20CR extends to 1851 in V2c,
		so that functionality is needed when handling that reanalysis data set. '''

		if repr_interval < 0:
			raise ValueError('Interval must be a positive integer or zero.')
		
		datestr = '%04d' % date.year
		if repr_interval >= 1:
			datestr += '%02d' % date.month
			if repr_interval >= 2:
				datestr += '%02d' % date.day
				if repr_interval >= 3:
					datestr += '%02d' % date.hour

		return datestr

	# Special cases: Only one datetime object given
	if type(dates) == dt:
		return date_format(dates, 3)
	
	elif len(dates) == 1:
		return date_format(dates[0], 3)
		
	# General case: Several datetime objects given
	dta = min(dates)
	dtz = max(dates)
	
	if dta.hour == 0:
		if dta.day == 1:
			if dta.month == 1:
				# year only
				repr_interval = 0
			else:
				# year and month
				repr_interval = 1
		else:
			# year, month and day
			repr_interval = 2
	else:
		# year, month, day and hour
		repr_interval = 3
	
	# Try to find an applicable time aggregation interval, either by argument or from the data set
	if agg:
		agg = getattr(tagg, agg)
		timestep = agg.interval
	elif se.conf.timestep:
		if type(se.conf.timestep) == str: 
			raise NotImplementedError()
		else:
			timestep = se.conf.timestep
			agg = tagg.get_by_interval(se.conf.timestep)

	else:
		agg = None

	dates = set(dates)
	
	dates = set(dates)
	
	# If an aggregation interval is available, find out if the given dates are contiguous 
	# and simplify the date representation accordingly
	if agg:
		agg = agg(dta, dtend=dtz+timestep, dtd=timestep)

		sep = '-'
		contiguous = True

		aggdates = list(agg)
		if not len(aggdates) == len(dates):
			sep = '..'
			contiguous = False
		else:
			for date in aggdates:
				if not date in dates:
					sep = '..'
					contiguous = False
					break
		
		# First time step of the next year 
		# -> First time step of the agg interval starting in the next year 
		# -> first time step of the last agg interval starting in the old year
		yearly = tagg.cal_year(dta)
		monthly = tagg.cal_month(dta)
		last4year = agg.start(yearly.start_next(aggdates[-1])) - timestep
		last4month = agg.start(monthly.start_next(aggdates[-1])) - timestep
		#
		if dtz.hour == last4month.hour:
			if dtz.day == last4month.day:
				if dtz.month == last4year.month:
					repr_interval = max(repr_interval, 0)
				else:
					repr_interval = max(repr_interval, 1)
			else:
				repr_interval = max(repr_interval, 2)
		else:
			repr_interval = 3

		fst = date_format(dta, repr_interval)
		lst = date_format(dtz, repr_interval)
		if fst == lst and contiguous:
			ret = fst
		else:
			ret = fst + sep + lst

	else:
		ret = date_format(dta, 3) + '..' + date_format(dtz, 3)

	return ret


def str2dts(periodstr, agg=None, epoch=None):
	''' Inverse of dts2str: Convert string to list of datetime objects '''

	known_formats = {'cal_year':'%Y', 'cal_month': '%Y%m', 'daily': '%Y%m%d', 'hourly': '%Y%m%d%H'}


	# Not defined in in function signature as they should be evaluated on run-time and not in import-time
	if not agg:
		agg = se.conf.timestep
	if not epoch:
		epoch = se.conf.epoch

	if periodstr.find('..') >= 0:
		raise ValueError('Cannot create list of dates for a period string representing a non-contiguous set of dates')

	if periodstr.find('-') >= 0:
		sep = periodstr.find('-')
		fst = periodstr[:sep]
		lst = periodstr[sep+1:]
	else:
		fst = periodstr
		lst = None
	
	for pagg, fmt in known_formats.items():
		try:
			fstp = dt.strptime(fst, fmt)
		except ValueError:
			continue
		else:
			break
	if not fstp:
		raise ValueError('The string `%s` could not be parsed into a time.' % (fst))
	
	agg = tagg.get_by_interval(agg)
	pagg = getattr(tagg, pagg)

	if lst:
		lstp = dt.strptime(lst, fmt)
	else:
		lstp = fstp
	
	lstp = pagg(epoch).start_next(lstp)
	
	dates = list(agg(fstp, lstp))

	return dates

# the end
