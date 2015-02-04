#!/usr/bin/env python

import os
import sys

import math
from datetime import datetime as dt, timedelta as td
import calendar
import numpy as np
from scipy.special import erfinv

import imp
pth = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
dynlib = imp.load_dynamic('dynlib', '%s/dynlib.so' % pth)

#
# Automatic scaling according to netcdf attributes "scale_factor" and "add_offset"
def scale(var, cut=(slice(None),slice(None),slice(None)), bench=False):
	if hasattr(var, 'scale_factor') or hasattr(var, 'add_offset'):
		if bench:
			begin = dt.now()
		# Python/numpy version is faster than Fortran function
		var_dat = var[cut]*getattr(var, 'scale_factor', 1.0) + getattr(var, 'add_offset', 0.0)
		#u_dat = dynlib.conv.scaleoff(u_dat, getattr(u, 'scale_factor', 1.0), getattr(u, 'add_offset', 0.0))
		if bench:
			print 'Python scaleoff', dt.now()-begin
	
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
	# Some map projections need the lon, lat array to be C-aligned.
	lon = np.ascontiguousarray(concat1(x))
	lat = np.ascontiguousarray(concat1(y))

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
			begin = dt.now()
		res = func(*args) 
		if bench:
			print 'Calculation', dt.now()-begin
	
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
# Rotate an angle from relative to the x-axis to relative to the wind direction
def rotate_natural(angle, u, v):
	angle -= np.arctan2(v,u)
	angle[angle < -np.pi/2] += np.pi
	angle[angle >= np.pi/2] -= np.pi

	return angle


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
				# python starts counting at zero, unlike fortran
				j = round(fronts[t,f,n,1] -1)
				i = round(fronts[t,f,n,0] -1) % s[1]
				masks[f][t,j,i] = True

	return masks

#
# return a 3d boolean array where line points are True, elsewere False
# OBS: Obsolete: use the Fortran version dynlib.utils.mask_lines instead
def mask_lines(lines, loff, s=(361,720)):
	mask = np.zeros((lines.shape[0], s[0], s[1]), dtype='bool')

	for t in range(lines.shape[0]):
		for n in range(loff[t].max()):
			# python starts counting at zero, unlike fortran
			j = round(lines[t,n,1] -1)
			i = round(lines[t,n,0] -1) % s[1]
			mask[t,j,i] = True

	return mask

#
# return a 3d boolean array where line points contain values, elsewere 0
def mask_lines_saveinfo(lines, loff, dat=None, s=(361,720)):
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

#
# return a 3d boolean array where line points are smoothly masked
def smear_lines(lines, loff, s=(361,720)):
	filtr_len = 5
	filtr_func = mk_gauss(0, 1)
	filtr = np.array(map(filtr_func, range(-filtr_len,filtr_len+1)))
	filtr /= sum(filtr)

	mask = dynlib.utils.mask_lines(s[1], s[0], lines, loff)
		
	return dynlib.utils.filter_xy(mask, filtr)



#
# mask parts of the field that are insignificant at the chosen level
def mask_insignificant(dat, mean, sig, nsig):
	dat[np.logical_and(dat < mean + nsig*sig, dat > mean - nsig*sig)] = np.nan

	return dat


# 
# Inverse of the CDF of the Gaussian distribution
igauss = lambda p: np.sqrt(2)*erfinv(2*p-1.0)

#
# Return a function calculating a Gaussian with the given mean (x0) and standard deviation (stddev)
def mk_gauss(x0,stddev):
	return lambda x: np.exp(-0.5*(x-x0)**2/stddev**2)/(np.sqrt(2*np.pi)*stddev)

#
# Calculate the most frequent value from a given histogram and bins
def cal_mfv(hist, bins):
	s = hist.shape[1:]
	mfv = np.zeros(s)
	for j in range(s[0]):
		for i in range(s[1]):
			bi = hist[:,j,i].argmax()
			mfv[j,i] = (bins[bi+1]+bins[bi])/2.0
	return mfv

#
# Generate (interpolation) points for cross section
def sect_gen_points(coords, m, dxy):
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

#
# General temporal aggregation 
def aggregate(dates, dat, agg):
	# We assume the time series to be equally spaced in time
	dtd = dates[1] - dates[0]

	# Helper functions generating functions :-)
	def first_func_gen(tdays):
		epoch = dt(1979,1,1)
		return lambda date: epoch + td(((date - epoch).days / tdays)*tdays)
	def last_func_gen(tdays):
		epoch = dt(1979,1,1)
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
	elif agg == 'cal_season':
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
	
	# Aggregate by ISO calendar weeks (which are aligned with 1979-1-1, so '7d' would be identical)
	elif agg == 'cal_weekly' or agg == '7d' or agg == 'weekly':
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
	
	# Aggregate by 10-day periods, starting 1979-1-1
	elif agg == '10d':
		first_func = first_func_gen(10)
		last_func = last_func_gen(10)
		agg = 'func'
	# Aggregate by 5-day periods, starting 1979-1-1
	elif agg == '5d':
		first_func = first_func_gen(5)
		last_func = last_func_gen(5)
		agg = 'func'
	# Aggregate by 3-day periods, starting 1979-1-1
	elif agg == '3d':
		first_func = first_func_gen(3)
		last_func = last_func_gen(3)
		agg = 'func'
	# Aggregate by 2-day periods, starting 1979-1-1
	elif agg == '2d':
		first_func = first_func_gen(2)
		last_func = last_func_gen(2)
		agg = 'func'
	# Aggregate by 24h periods, starting 1979-1-1
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
# Time series smoothing (running mean)
def t_smooth(ts, smooth):
	sts = np.zeros(ts.shape)
	
	a = -(smooth-1)/2
	z = (smooth-1)/2
	sts[:z] = np.nan
	sts[a:] = np.nan
	for i in range(a, z+1):
		tsslc = slice(z+i,a+i)
		if a+i == 0:
			tsslc = slice(z+i,None)
		else:
			tsslc = slice(z+i,a+i)

		sts[z:a] += ts[tsslc]
	
	sts /= smooth

	return sts


#
# Varimax rotation for EOFs as introduced in Kaiser (1958)
#  -> if raw = False use normal varimax, otherwise raw varimax
def varimax(phi, raw=False, gamma = 1.0, q = 20, tol = 1e-6):
	if not raw:
		norm = np.sqrt((phi**2).sum(axis=1))
		phi /= norm[:,np.newaxis]

	p, k = phi.shape
	R = np.eye(k)
	d = 0
	for i in xrange(q):
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
