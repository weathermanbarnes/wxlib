#!/usr/bin/env python
# -*- encoding: utf-8

import os
import sys

import numpy as np
import scipy.io.netcdf as nc
import scipy.io.matlab as mat
import pytz

from settings import conf as c
import utils
from gridlib import grid_by_nc, grid_by_static

from datetime import datetime as dt, timedelta as td

import imp
pth = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
dynlib = imp.load_dynamic('dynlib', '%s/dynlib.so' % pth)
dynlib_version = (''.join(dynlib.consts.version)).strip()

# #############################################################################
# 1. Open all sorts of files
# 

# Find and open files by filename (without ending!)
def metopen(filename, q, cut=slice(None), verbose=False, no_dtype_conversion=False, no_static=False):
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
			dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.npz'
		elif os.path.exists(path+'/'+filename+'.mat'):
			f   = mat.loadmat(path+'/'+filename+'.mat')
			dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.mat'
		elif os.path.exists(path+'/'+filename+'.nc'):
			f   = nc.netcdf_file(path+'/'+filename+'.nc', 'r')
			var = f.variables[q]
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
	s = dat.shape
	if not len(s) == 3 or not s[1:] == (361,720):
		raise NotImplementedError, 'dat does not seem to be a ERA-Interim like (t,y,x)-array.'
	
	now = dt.now(pytz.timezone('Europe/Oslo'))
	of = nc.netcdf_file(c.opath+'/'+(c.file_std % {'time': time, 'plev': plev, 'q': c.qi[q]})+'.nc', 'w')
	of._attributes = {'Conventions': 'CF-1.0', 
			'history': '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)
	}

	of.createDimension('time', dat.shape[0])
	of.createDimension('latitude', dat.shape[1])
	of.createDimension('longitude', dat.shape[2])

	ot = of.createVariable('time', 'i', ('time',))
	ot._attributes = {'long_name': 'time', 'units': static.t_unit}
	ot[:] = static.t
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
	llen = int(datoff.max())
	# The maxmimum amount of lines for all time steps
	olen = int(np.max([datoff_t.argmin() for datoff_t in datoff[:,1:]]))

	of.createDimension('time', dat.shape[0])
	of.createDimension('pointindex', llen)
	of.createDimension('infotype', 3)
	of.createDimension('lineindex', olen)

	ot = of.createVariable('time', 'i', ('time',))
	ot._attributes = {'long_name': 'time', 'units': static.t_unit}
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
	fo, oro = metopen('static', 'oro', cuts, verbose, no_dtype_conversion, True)
	static = grid_by_static(fo, cut=cuts)
	static.oro = oro[::]
	fo.close()

	return static


# #############################################################################
# 2. Generalised data fetchers
# 

# Generalised data fetcher for instantaneous or short-term averaged fields
def get_instantaneous(q, dates, plevs=None, yidx=None, xidx=None, tavg=False, quiet=False, force=False, **kwargs):
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
				f, d, static = metopen(c.file_std % {'time': year, 'plev': plev, 'q': c.qi[q]}, q, cut=cut, **kwargs)
				s = tuple([1+tsmax-tsmin,len(plevs)] + list(d.shape)[1:])
				dat = np.empty(s, dtype=d.dtype)
			else:
				f, d = metopen(c.file_std % {'time': year, 'plev': plev, 'q': c.qi[q]}, q, cut=cut, no_static=True, **kwargs)
			dat[datcut,i,::] = d
			i += 1

	# Time-averaging if specified
	if tavg and len(dates) > 1:
		dat = dat.mean(axis=0)
	
	dat = dat.squeeze()
	
	return dat, static


# Get aggregated (average, standard deviation, etc.) fields
def get_aggregate(q, year=None, plev=None, yidx=None, xidx=None):
	raise NotImplementedError, 'If you need it, implement it!'

	return dat


# the end
