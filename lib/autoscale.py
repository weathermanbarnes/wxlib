#!/usr/bin/env python
# -*- encoding: utf-8


''' A more intelligent and configurable autoscaling

The standard autoscaling by matplotlib has several deficiencies.

 #. It is based on the minimum and maximum of the field. As a result, spiky
    data will screw up the scaling completely. 
 #. It does not take the extend-keyword into account.
 #. It can result in very unintuitive intervals like 1.25 to whatever power is appropriate.
 #. It cannot enforce scales that are symmetric around zero.

This module provides an alternative mechanism to automatically determine
appropriate scales that alleviates all these deficiencies. The autoscaling
is configurable through several configuration keys in the plotconfig.
'''


import numpy as np


DEFAULT_STEPS = 7
DEFAULT_INTERVALS = np.array([1,2,3,5,10])
DEFAULT_PERCENTILES = (0.02, 0.98)

def _find_scale(lthres, uthres, target_steps, symmetric_zero, intervals=DEFAULT_INTERVALS, intervals_periodic=True):
	''' Find the best matching scale
	
	It takes into account the given intervals and the given target amount of steps. 
	'''

	if not type(intervals) in [np.ndarray, list, tuple]:
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


def _find_thres(dat, mask, extend, exceed_percentiles, symmetric_zero):
	''' Find upper and lower limit of the colorbar
	
	The Result is based on the extend keyword and given exceedence percentiles. 
	'''

	nanmask = np.logical_not(np.isnan(dat[mask]))
	dat = sorted(dat[mask][nanmask].flatten())
	ldat = len(dat)

	if len(dat) == 0:
		return -0.5, 0.5
	if min(dat) == max(dat):
		return min(dat)-0.5, max(dat)+0.5

	lower, upper = exceed_percentiles

	if extend in ['both', 'lower', 'min']:
		lthres = dat[int(np.floor(lower*ldat))]
	else:
		lthres = dat[0] 
	if extend in ['both', 'upper', 'max']:
		uthres = dat[int(np.floor(upper*ldat))]
	else:
		uthres = dat[-1]
	
	if symmetric_zero:
		lthres = -1*max(abs(lthres), abs(uthres))
		uthres =    max(abs(lthres), abs(uthres))
	
	return lthres, uthres


def autoscale(dat, **conf):
	''' Determine the scale based on the data and config

	Parameters
	----------
	dat : np.ndarray
	    Data determining the scale
	
	Keyword arguments
	-----------------
	mask : np.ndarray with dtype bool
	    Default: No masking. Allows to take only parts of the data array into account. 
	extend : str
	    **Required**. Same as the matplotlib keyword argument ``extend``. Valid values are 
	    ``'neither'``, ``'min'``, ``'max'`` and  ``'both'``. 
	scale_exceed_percentiles : 2-tuple of float
	    Default: ``%s``. If ``extend`` allows data outside the main colorbar range, how 
	    much of the data is allowed to lie outside the main range? The first value is the 
	    maxmimum percentile for the first tick at the colorbar, the second value the minimum
	    percentile for the last tick at the colorbar.
	scale_symmetric_zero : bool
	    Default: ``False``. Force the scale to be symmetric around zero.
	scale_target_steps : int
	    Default: ``%d``. Number of colors in the main colorbar range. **Note**: The actual number 
	    might deviate slightly from the target number given here, if required by the data range 
	    and the allowed intervals.
	scale_intervals : np.ndarray
	    Default: ``%s``. Which tick intervals are allowed? 
	scale_interval_periodic : bool
	    Default: ``True``. If true, the numbers given in ``scale_intervals`` taken independently 
	    of the exponent, so an entry ``5`` allows a tick interval of ``0.5`` as well 
	    as one of ``5000``.
	
	Returns
	-------
	np.ndarray
	    Scale
	''' 
	
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


# Insert default values into docstring
autoscale.__doc__ %= (DEFAULT_PERCENTILES, DEFAULT_STEPS, str(DEFAULT_INTERVALS))


# the end
