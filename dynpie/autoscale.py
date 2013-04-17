#!/usr/bin/env python
# -*- encoding: utf-8
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  		DynLib -- A more intelligent and configurable autoscaling
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
import numpy as np

DEFAULT_STEPS = 7
DEFAULT_INTERVALS = np.array([1,2,3,5,10])
DEFAULT_PERCENTILES = (0.02, 0.98)

# Find the best matching scale, taking into account the given choice intervals and the given target amount of steps.
def _find_scale(lthres, uthres, target_steps, symmetric_zero, intervals=DEFAULT_INTERVALS, intervals_periodic=True):
	if not type(intervals) == np.ndarray and not type(intervals) == list:
		interval = intervals
	elif len(intervals) == 1:
		interval = intervals[0]
	elif target_steps:
		log10 = np.log(10)
		log10intervals = np.log(intervals)/log10
		log10interval  = np.log((uthres - lthres)/target_steps)/log10
		
		if intervals_periodic:
			exp = np.floor(log10interval)
			log10interval = log10interval % 1.0
		else:
			exp = 0
		
		i = 0
		diff = 9999.0
		while i < len(log10intervals) and diff > abs(log10intervals[i]-log10interval):
			diff = abs(log10intervals[i]-log10interval)
			i += 1
		
		i -= 1
		interval = intervals[i] * 10**int(exp)
	else:
		raise ValueError, 'Either interval or target_steps must be given'

	if symmetric_zero:
		lthres = (np.floor(lthres/interval - 0.5)+0.5)*interval
		uthres = (np.ceil(uthres/interval + 0.5)-0.5)*interval + interval/10.0
	else:
		lthres =  np.floor(lthres/interval)*interval
		uthres =  np.ceil(uthres/interval)*interval + interval/10.0

	return np.arange(lthres,uthres,interval)


# Find upper and lower limit of the colorbar, based on the extend keyword and given exceedence percentiles
def _find_thres(dat, mask, extend, exceed_percentiles, symmetric_zero):
	dat = sorted(dat[mask].flatten())
	ldat = len(dat)

	lower, upper = exceed_percentiles

	if extend == 'both' or extend == 'lower':
		lthres = dat[int(np.floor(lower*ldat))]
	else:
		lthres = dat[0] 
	if extend == 'both' or extend == 'upper':
		uthres = dat[int(np.floor(upper*ldat))]
	else:
		uthres = dat[-1]
	
	if symmetric_zero:
		lthres = -1*max(abs(lthres), abs(uthres))
		uthres =    max(abs(lthres), abs(uthres))
	
	return lthres, uthres


# Main function: takes the data and config; return the scale
def autoscale(dat, **conf):
	mask   = conf.get('mask', slice(None))
	extend = conf.get('extend')
	exceed = conf.get('scale_exceed_percentiles', DEFAULT_PERCENTILES)
	symzer = conf.get('scale_symmetric_zero', False)
	steps  = conf.get('scale_target_steps', DEFAULT_STEPS)
	intvls = conf.get('scale_intervals', DEFAULT_INTERVALS)
	intper = conf.get('scale_intervals_periodic', True)

	lthres, uthres = _find_thres(dat, mask, extend, exceed, symzer)
	scale = _find_scale(lthres, uthres, steps, symzer, intvls, intper)

	return scale


# the end
