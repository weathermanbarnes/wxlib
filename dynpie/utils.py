#!/usr/bin/env python

import math
import datetime
import numpy as np
from scipy.special import erfinv

#
# Automatic scaling according to netcdf attributes "scale_factor" and "add_offset"
def scale(var, cut=(slice(None),slice(None),slice(None)), bench=False):
	if hasattr(var, 'scale_factor') or hasattr(var, 'add_offset'):
		if bench:
			begin = datetime.datetime.now()
		# Python/numpy version is faster than Fortran function
		var_dat = var[cut]*getattr(var, 'scale_factor', 1.0) + getattr(var, 'add_offset', 0.0)
		#u_dat = dynlib.conv.scaleoff(u_dat, getattr(u, 'scale_factor', 1.0), getattr(u, 'add_offset', 0.0))
		if bench:
			print 'Python scaleoff', datetime.datetime.now()-begin
	
	else:
		var_dat = var[cut]
	
	return var_dat

#
# Reduce to i2 values with add_offset and scale_factor
def unscale(var):
	maxv  = var.max()
	minv  = var.min()
	# divide in 2^16-2 intervals, values from -32766 -> 32767 ; reserve -32767 as missing value
	scale = (maxv-minv)/65534.0
	off   = +32766.5*scale + minv

	res = np.round((var[::] - off)/scale)

	return res.astype('>i2'), scale, off

#
# Concatenate one latitude band in x-direction, taking over the values of 
# the first latitude band to emulate a cyclic field in Basemap plots
def concat1(data):
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

#
# Concatenate one latitude band in x-direction, taking over the values of 
# the first latitude band to emulate a cyclic field in Basemap plots
def concat1lonlat(x, y):
	# Some map projections need the lon, lat array to be Fortran-aligned.
	lon = np.asfortranarray(concat1(x))
	lat = np.asfortranarray(concat1(y))

	lon[:,-1] += 360.0

	return lon, lat

#
# Generic calculation preparations and the actual call of the Fortran function
def call(func, vars, grid, cut=(slice(None),slice(None),slice(None)), bench=False):
	if not grid.nz or grid.nz == 1:
		print '3D mode'
		if len(vars[0].shape) == 4:
			raise NotImplementedError
		args = []
		for var in vars:
			if len(var.shape) == 3:
				if getattr(var, 'dtype', None) in ('i2','>i2','<i2'):
					args.append(scale(var, cut, bench=bench))
				else:
					args.append(var[cut])
			else:
				args.append(var) 
			
		args.extend([grid.dx[cut[1:]], grid.dy[cut[1:]]])
		if bench:
			begin = datetime.datetime.now()
		res = func(*args) 
		if bench:
			print 'Calculation', datetime.datetime.now()-begin
	
	elif not grid.nt or grid.nt == 1:
		print '2D mode'
		raise NotImplementedError
	
	else:
		print '4D mode'
		raise NotImplementedError
		#for t in len(u.shape[0]):
		#	dylib.diag.def(u[t,:,:,:], v[t,:,:,:], grid.dx, grid.dy)
	
	return res


#
# Reimplementation of the recpective function in dynlib.diag for benchmarking.
def def_angle(u_dat, v_dat, grid):
	deff = np.zeros(u_dat.shape)
	for k in range(u_dat.shape[0]):
		for j in range(1,grid.ny-1):
			for i in range(1,grid.nx-1):
				def_shear   = (u[k,j+1,i]-u[k,j-1,i])/grid.dy[j,i] \
					    + (v[k,j,i+1]-v[k,j,i-1])/grid.dx[j,i]
				def_stretch = (u[k,j,i+1]-u[k,j,i-1])/grid.dx[j,i] \
					    - (v[k,j+1,i]-v[k,j-1,i])/grid.dy[j,i]
				deff[k,j,i] = 0.5*math.atan2(def_shear, def_stretch)
			if grid.cyclic_ew:
				def_shear   = (u[k,j+1,i]-u[k,j-1,i])/grid.dy[j,1] \
					    + (v[k,j,  1]-u[k,j, -1])/grid.dx[j,1]
				def_stretch = (u[k,j,  1]-u[k,j, -1])/grid.dx[j,1] \
					    - (v[k,j+1,i]-v[k,j-1,i])/grid.dy[j,1]
				deff[k,j,1] = 0.5*math.atan2(def_shear, def_stretch)
				def_shear   = (u[k,j+1,i]-u[k,j-1,i])/grid.dy[j,1] \
					    + (v[k,j,  0]-u[k,j, -2])/grid.dx[j,1]
				def_stretch = (u[k,j,  0]-u[k,j, -2])/grid.dx[j,1] \
					    - (v[k,j+1,i]-v[k,j-1,i])/grid.dy[j,1]
				deff[k,j,-1] = 0.5*math.atan2(def_shear, def_stretch)

	return deff


#
# Unflatten a flattened front array, using the froff list; 
#  separately for cold/warm and stationary fronts
def unflatten_fronts(fronts, froff, minlength=1):
	cold = [__unflatten_fronts_t(fronts[t,0,:,:], froff[t,0,:], minlength) for t in range(fronts.shape[0])]
	warm = [__unflatten_fronts_t(fronts[t,1,:,:], froff[t,1,:], minlength) for t in range(fronts.shape[0])]
	stat = [__unflatten_fronts_t(fronts[t,2,:,:], froff[t,2,:], minlength) for t in range(fronts.shape[0])]

	return cold, warm, stat

# unflatten one time step and front type
def __unflatten_fronts_t(fronts, froff, minlength):
	fronts = [fronts[froff[i]:froff[i+1],:] for i in range(froff.shape[0]-1)]
	fronts = filter(lambda lst: len(lst) >= minlength, fronts)

	return fronts

#
# return a 3d boolean array where frontal points are True, elsewere False
def mask_fronts(fronts, froff, s=(361,720)):
	masks = [np.zeros((len(fronts), s[0], s[1]), dtype='bool'),
		 np.zeros((len(fronts), s[0], s[1]), dtype='bool'),
		 np.zeros((len(fronts), s[0], s[1]), dtype='bool') ]

	for t in range(len(fronts)):
		for f in range(3):
			for n in range(froff[t,f].max()):
				j = round(fronts[t,f,n,1])
				i = round(fronts[t,f,n,0])
				masks[f][t,j,i] = True

	return masks


# 
# Inverse of the CDF of the Gaussian distribution
igauss = lambda p: np.sqrt(2)*erfinv(2*p-1.0)

#
# Calculate the most frequent value from a given histogram and bins
def cal_mfv(hist, bins):
	mfv = np.zeros(s)
	for j in range(s[0]):
		for i in range(s[1]):
			bi = hist[:,j,i].argmax()
			mfv[j,i] = (bins[bi+1]+bins[bi])/2.0
	return mfv
#
