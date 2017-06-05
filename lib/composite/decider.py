#!/usr/bin/env python
# -*- encoding: utf-8

''' This module defines decider functions for composites

Each decider implements a criterion to deterime if a given time step should 
be part of the composite or not. The decision can be based on either time
itself, a provided time series or test data.

Several deciders can be combined by logical operators. For example a combination
of a decider for positive NAO phases and another decider for winter can be used
to construct a combined decider for positive NAO during winter.

Furthermore, the module provides functions to create lagged composites from a
given decider, and to create all combinations between two lists of composites
(e.g. for combining all variability indexes with all seasons).
'''

from __future__ import absolute_import, unicode_literals, print_function

from copy import deepcopy

from ..settings import conf
from ..shorthands import np, dt, td, metopen


# General base class
class decider(object):
	''' Abstract base class for all deciders
	
	The base class defines the shared API for deciders implements the logical 
	operators to combine deciders. 

	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	q : str
	    *Optional*. If applicable, the variable name identifier for required test 
	    data, following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``. 
	plev : str
	    *Optional*. If applicable, the vertical level of the required test data, 
	    following the ECMWF conventions, e.g. ``'700'`` for 700 hPa or ``'pv2000'`` 
	    for the PV2-surface.
	'''

	def __init__(self, name, q=None, plev=None, rotation_center=None):
		self.name = name
		self.rotation_center = rotation_center
		
		if q and plev:
			self.required_qs = (plev, q)

		if plev:
			self.plev = plev
			self.name += '_%s' % plev
		if q:
			self.q = q
			self.name += '_%s' % q

		return

	
	def __combine(self, b, ret):
		''' Combines two decider, incorporates anything shared between __or__ and __and__ '''

		ret.reset = lambda : self.reset() | b.reset()
		
		selfrot = not type(self.rotation_center) == type(None)
		brot = not type(b.rotation_center) == type(None)

		if selfrot and brot:
			if not self.rotation_center == b.rotation_center:
				raise TypeError('Cannot combine two deciders with different rotation centers')
			ret.rotation_center = self.rotation_center

		elif selfrot:
			ret.rotation_center = self.rotation_center

		elif brot:
			ret.rotation_center = b.rotation_center

		return ret

	def __or__(self, b):
		''' Implements the ``deciderA | deciderB`` syntax '''

		ret = decider(self.name+b.name)
		ret.match = lambda date, tidx, prv, cur, nxt: self.match(date, tidx, prv, cur, nxt) | b.match(date, tidx, prv, cur, nxt)

		ret = self.__combine(b, ret)

		return ret

	def __and__(self, b):
		''' Implements the ``deciderA & deciderB`` syntax '''

		ret = decider(self.name+'@'+b.name)
		ret.match = lambda date, tidx, prv, cur, nxt: self.match(date, tidx, prv, cur, nxt) & b.match(date, tidx, prv, cur, nxt)

		ret = self.__combine(b, ret)

		return ret
	
	# The deciding function, returns True/False
	# to be overriden by derived classes.
	def match(self, date, tidx, prv, cur, nxt):
		''' The actual deciding function if a given time should be part of the composite

		Parameters
		----------
		date : datetime
		    Date of the given time
		tidx : int
		    Time index of the given time
		prv : dict of np.ndarray. Keys is the dict are 2-tuples with the entries
		    (variable identifer, vertical level).
		    Test data for the previous year. Keys is the dict are 2-tuples with the entries
		    (variable identifer, vertical level).
		cur : dict of np.ndarray
		    Test data for the current year. Keys is the dict are 2-tuples with the entries
		    (variable identifer, vertical level).
		next : dict of np.ndarray
		    Test data for the next year. Keys is the dict are 2-tuples with the entries
		    (variable identifer, vertical level).
	
		Returns 
		-------
		bool
		    True of the given time should be part of the composite, and False otherwise.
		'''

		return False
	
	def reset(self):
		''' Reset the decider for a new pass through time 
		
		Returns
		-------
		False
		    Important for resetting combined deciders.
		'''

		return False


# Criterion: exceedance of a lower threshold at a given location
class dat_lowerbound(decider):
	''' Decider based on the positive exceedence of a threshold, applied to test data 
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	q : str
	    *Optional*. If applicable, the variable name identifier for required test 
	    data, following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``. 
	plev : str
	    *Optional*. If applicable, the vertical level of the required test data, 
	    following the ECMWF conventions, e.g. ``'700'`` for 700 hPa or ``'pv2000'`` 
	    for the PV2-surface.
	pos : 2-tuple of int or 2-tuple of slice
	    Coordinates of the point or slice to test, in grid point indexes.
	thres : float
	    Threshold.
	agg : callable
	    *Optional*, default: median. A function aggregating several values into one.
	    Relevant if pos contained slices to aggregate the extraced grid area into one
	    value to be tested against the threshold.
	'''

	def __init__(self, name, q, plev, pos, thres, agg=np.median):
		decider.__init__(self, name, q, plev)
		self.thres = thres
		self.yidx, self.xidx = pos
		self.agg = agg

		return

	def _get_vals(self, tidx, prv, cur, nxt):
		''' Extract the required data from the diven test data dict '''

		if tidx < 0:
			if (self.plev, self.q) not in prv:
				return None
			vals = prv[self.plev, self.q][tidx, self.yidx, self.xidx]
		elif tidx >= cur[self.plev, self.q].shape[0]:
			if (self.plev, self.q) not in nxt:
				return None
			vals = nxt[self.plev, self.q][tidx - cur[self.plev, self.q].shape[0], self.yidx, self.xidx]
		else:
			vals = cur[self.plev, self.q][tidx,self.yidx,self.xidx]

		return vals

	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = decider.match.__doc__

		vals = self._get_vals(tidx, prv, cur, nxt)
		if type(vals) == type(None):
			return False

		return self.agg(vals) >= self.thres
	

# Criterion: (negative) exceedance of an upper threshold at a given location
class dat_upperbound(dat_lowerbound):
	''' Decider based on the negative exceedence of a threshold, applied to test data 
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	q : str
	    *Optional*. If applicable, the variable name identifier for required test 
	    data, following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``. 
	plev : str
	    *Optional*. If applicable, the vertical level of the required test data, 
	    following the ECMWF conventions, e.g. ``'700'`` for 700 hPa or ``'pv2000'`` 
	    for the PV2-surface.
	pos : 2-tuple of int or 2-tuple of slice
	    Coordinates of the point or slice to test, in grid point indexes.
	thres : float
	    Threshold.
	agg : callable
	    *Optional*, default: median. A function aggregating several values into one.
	    Relevant if pos contained slices to aggregate the extraced grid area into one
	    value to be tested against the threshold.
	'''

	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = dat_lowerbound.match.__doc__

		vals = self._get_vals(tidx, prv, cur, nxt)
		if type(vals) == type(None):
			return False

		return self.agg(vals) < self.thres


# Criterion: boolean array true at a given location
class dat_boolean(dat_lowerbound):
	''' Decider based on the negative exceedence of a threshold, applied to test data 
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	q : str
	    *Optional*. If applicable, the variable name identifier for required test 
	    data, following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``. The test
	    data must be boolean variable.
	plev : str
	    *Optional*. If applicable, the vertical level of the required test data, 
	    following the ECMWF conventions, e.g. ``'700'`` for 700 hPa or ``'pv2000'`` 
	    for the PV2-surface.
	pos : 2-tuple of int or 2-tuple of slice
	    Coordinates of the point or slice to test, in grid point indexes.
	thres : float
	    Threshold.
	agg : callable
	    *Optional*, default: any. A function aggregating several boolean values into 
	    one. Relevant if pos contained slices to aggregate the extraced grid area into 
	    one boolean value to be returned.
	'''

	def __init__(self, name, q, plev, pos, agg=np.any):
		decider.__init__(self, name, q, plev)
		self.yidx, self.xidx = pos
		self.agg = agg

		return
	
	# Decider
	def match(self, date, tidx, prv, cur, nxt):
		vals = self._get_vals(tidx, prv, cur, nxt)
		if type(vals) == type(None):
			return False

		return self.agg(vals)
	

# Criterion: exceedance of a lower threshold in a given time series
class ts_lowerbound(decider):
	''' Decider based on the positive exceedence of a threshold, applied to a time series
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	ts : dict
	    Time series dictionary, containing a list of datetime objects and a list of values.
	thres : float
	    Threshold.
	'''
	
	def __init__(self, name, ts, thres):
		decider.__init__(self, name)
		self.tidx   = 0
		self.dates  = ts['dates']
		self.values = ts['values']
		self.thres  = thres

		return
	
	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = decider.match.__doc__

		if date < self.dates[0] or date >= self.dates[-1]:
			return False

		while self.dates[self.tidx+1] <= date:
			self.tidx += 1
		
		return self.values[self.tidx] >= self.thres

	def reset(self):
		__doc__ = decider.reset.__doc__

		self.tidx = 0

		return False


# Criterion: (negative) exceedance of an upper threshold in a given time series
class ts_upperbound(ts_lowerbound):
	''' Decider based on the negative exceedence of a threshold, applied to a time series
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	ts : dict
	    Time series dictionary, containing a list of datetime objects and a list of values.
	thres : float
	    Threshold.
	'''

	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = ts_lowerbound.match.__doc__

		if date < self.dates[0] or date >= self.dates[-1]:
			return False

		while self.dates[self.tidx+1] <= date:
			self.tidx += 1
		
		return self.values[self.tidx] < self.thres


# Criterion: Given time series equal to given value
class ts_equal(ts_lowerbound):
	''' Decider checking a time series is equal to a given value

	The time series should consist of integers to avoid rounding problems.
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	ts : dict
	    Time series dictionary, containing a list of datetime objects and a list of values.
	thres : float
	    Given value to be tested if the time series is equal to.
	'''

	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = ts_lowerbound.match.__doc__

		if date < self.dates[0] or date >= self.dates[-1]:
			return False

		while self.dates[self.tidx+1] <= date:
			self.tidx += 1
		
		return self.values[self.tidx] == self.thres


# Criterion: Given time series equal to given value
class ts_notnan(ts_lowerbound):
	''' Decider checking if a time series is not NaN
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	ts : dict
	    Time series dictionary, containing a list of datetime objects and a list of values.
	'''
	
	def __init__(self, name, ts):
		ts_lowerbound.__init__(self, name, ts, None)

		return
	
	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = ts_lowerbound.match.__doc__

		if date < self.dates[0] or date >= self.dates[-1]:
			return False

		while self.dates[self.tidx+1] <= date:
			self.tidx += 1

		return not np.isnan(self.values[self.tidx])


# Criterion: Match a given month
class month(decider):
	''' Decider checking if in a given month 
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	month : int
	    Month of the year.
	'''

	def __init__(self, name, month):
		decider.__init__(self, name)
		self.month = month

		return
	
	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = decider.match.__doc__

		return date.month == self.month


# Criterion: Match a given time of the day
class timeofday(decider):
	''' Decider checking if within a given range of hours 
	
	Parameters
	----------
	name : str
	    Name of the composite, to be used for saving.
	begin_hour : int
	    Start hour of the time interval, the hour sharp is included.
	end_hour : int
	    End hour of the time interval, the hour sharp is excluded.
	'''
	
	def __init__(self, name, begin_hour, end_hour):
		decider.__init__(self, name)
		self.begin = begin_hour
		self.end = end_hour

		return
	
	# Decider
	def match(self, date, tidx, prv, cur, nxt):
		__doc__ = decider.match.__doc__

		return date.hour >= self.begin and date.hour < self.end


def __timelag_one(decider, lidx, dtidx, dt):
	''' Create a time-lagged version of a decider '''

	tl_decider = deepcopy(decider)
	tl_decider.__match = tl_decider.match
	def match(date, tidx, prv, cur, nxt):
		return tl_decider.__match(date-dt, tidx-dtidx, prv, cur, nxt)
	tl_decider.match = match

	if dtidx >= 0: sgn = '+'
	else: sgn = '-'
	tl_decider.name += '_i%02dd%s%02d' % (lidx, sgn, abs(dtidx))

	return tl_decider

def timelag(decider, dtidxs, tstep=conf.timestep):
	''' Create a series of time-lagged versions of one decider 

	Parameters
	----------
	decider : decider
	    Decider to be time-lagged.
	dtidxs : list of int
	    List of time steps by which the decider should be time-lagged.
	tstep : timedelta or relativedelta
	    *Optional*, default ``conf.timestep``. Time step of the data set.

	Returns
	-------
	dict of list of decider
	    List of time-lagged versions of the decider, wrapped in a named group
	    such that they will by default be saved together.
	'''

	tl_deciders = {decider.name: map(lambda lidx, dtidx: __timelag_one(decider, lidx, dtidx, dtidx*tstep), 
			zip(range(len(dtidxs)), dtidxs))}

	return tl_deciders

def matrix(list1, list2):
	''' Create a list of deciders containing all combinations of the given decider lists
	
	If the two lists are, for example, a list of climate variability indexes and a list 
	of seasons, a list of deciders for all variability indexes for all seasons is returned.

	Parameters
	----------
	list1: list of decider
	    First list of deciders to be combined.
	list2: list of decider
	    Second list of deciders to be combined. Each decider in this list makes a resulting category.

	Returns
	-------
	list of decider
	    List of combinations.
	'''

	# Flatten dicts
	if type(list1) == dict:
		list1 = [item for group in list1.values() for item in group]
	if type(list2) == dict:
		list2 = [item for group in list2.values() for item in group]
	if len(list1) == 0 or len(list2) == 0:
		return {}
	
	# Wrap single deciders in a list
	if not type(list1) == list: list1 = [list1, ]
	if not type(list2) == list: list2 = [list2, ]

	return {item2.name: [item1 & item2 for item1 in list1] for item2 in list2}

# the end
