#!/usr/bin/env python
# -*- encoding: utf-8

import os.path
import numpy as np
import scipy.io.netcdf as nc
import scipy.io.matlab as mat
from settings import conf as c
import utils
from gridlib import grid_by_nc, grid_by_static


# #############################################################################
# 1. Open all sorts of files
# 

# Find and open files by filename (without ending!)
def metopen(filename, q, cut=slice(None), verbose=False, no_dtype_conversion=False, no_static=False):
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
			static = grid_by_nc(f)
			# TODO: Where to search for topography in nc files?
			static.oro = np.zeros((static.ny, static.nx))
			print 'Found '+path+'/'+filename+'.nc'
		else:
			continue
		
		if not no_dtype_conversion:
			dat = dat.astype('f8')

		if not no_static:
			if not static:
				if type(cut) == tuple:
					if len(cut) > 1:
						cuts = tuple(list(cut)[1:])
					else: 
						cuts = slice(None)
				else: 
					cuts = slice(None)
				static = get_static(cuts, verbose, no_dtype_conversion)
						
		else:
			return f, dat

		return f, dat, static
	
	raise RuntimeError, '%s.* not found in any data location.' % filename


def get_static(cuts=slice(None), verbose=False, no_dtype_conversion=False):
	fo, oro = metopen('static', 'oro', cuts, verbose, no_dtype_conversion, True)
	static = grid_by_static(fo)
	static.oro = oro[::]
	fo.close()

	return static


# #############################################################################
# 2. Generalised data fetchers
# 

# Generalised data fetcher for instantaneous or short-term averaged fields
def get_instantaneous(q, dates, plevs=None, yidx=None, xidx=None, tavg=True, quiet=False):
	# None means "take everything there is as pressure levels"
	if not plevs:
		plevs = c.plevs
	else:
		plevs = [plevs,]
	
	if yidx == None:
		yidxs = slice(None)
	else:
		yidxs = yidx
	
	if xidx == None:
		xidxs = slice(None)
	else:
		xidxs = xidx

	# Convert dates to time indexes
	if type(dates) not in ([np.ndarray, list, tuple, set]):
		dates = [dates, ]
	tidxs = map(lambda x: (x.timetuple().tm_yday-1)*4 + int(x.hour/6), dates)

	# Construct the slice
	cut = (slice(min(tidxs),max(tidxs)+1), yidxs, xidxs)
	
	# One ore more vertical levels?
	if len(plevs) > 1:
		i = 0
		#dat = dat.squeeze()
		for plev in plevs:
			if not quiet:
				print "Reading from "+c.file_std % (dates[0].year, plev, q)
			if plev == plevs[0]:
				f, d, static = metopen(c.file_std % (dates[0].year, plev, q), c.q[q], cut=cut, no_static=True)
				dat = np.zeros((1+max(tidxs)-min(tidxs), len(c.plevs), d.shape[1], d.shape[2]))
				
			else:
				f, d = metopen(c.file_std % (dates[0].year, plev, q), c.q[q], cut=cut)
			dat[:,i,::] = d
			i += 1
	else:
		if not quiet:
			print "Reading from "+c.file_std % (dates[0].year, plevs[0], q)
		f, dat, static = metopen(c.file_std % (dates[0].year, plevs[0], q), c.q[q], cut=cut)
	
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
