#!/usr/bin/env python
# -*- encoding: utf-8

import os.path
import numpy as np
import scipy.io.netcdf as nc
import scipy.io.matlab as mat
from settings import conf as c
import utils
from gridlib import grid_by_nc, grid_by_static

from datetime import datetime as dt, timedelta as td


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
