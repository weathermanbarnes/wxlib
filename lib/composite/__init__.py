#!/usr/bin/python
# -*- encoding: utf-8

from .. import utils, metio
from ..shorthands import np, dt, td, metopen, metsave_timeless
from ..settings import conf

from scipy.io import savemat


# TODO: LINES should be centralised somewhere in the variable definitions!
LINES = {'fronts': 'froff', 'convl': 'cloff', 'defl': 'dloff', 'vorl': 'vloff', 'jetaxis': 'jaoff'}


cnt = ('cnt', None, 'Number of time steps contributing to the composite', '1')
conf.register_variable([cnt, ], [])


def _add(name, plev, q, dat, mean, hist):
	''' Add one time step to the composite '''

	if q in conf.q_bins:
		for bi in range(len(conf.q_bins[q])-1):
			upper = conf.q_bins[q][bi+1]
			lower = conf.q_bins[q][bi]
			if upper > lower:
				hist[name,plev,q][bi,(dat >= lower).__and__(dat <  upper)] += 1
			else:
				hist[name,plev,q][bi,(dat <  upper).__or__(dat >= lower)] += 1
	else:
		mean[name,plev,q][:,:] += dat

	return


def _get_testdat(time, test_qs, readhooks):
	''' Retrieve the data used in deciders for testing '''

	testdat = {}
	for test_plev, test_q in test_qs:
		if not test_q:
			continue
	
		if test_q in readhooks:
			f, testdat[test_plev, test_q] = readhooks[test_q](time, test_plev, test_q)
		else:
			f, testdat[test_plev, test_q] = metopen(conf.file_std % {'time': time, 'plev': test_plev, 'qf': conf.qf[test_q]}, test_q, no_static=True)

		if test_q in LINES:
			testdat[test_plev, test_q] = utils.mask_lines(testdat[test_plev, test_q], f.variables[LINES[test_q]][::])
	
	return testdat


def build(qs, tests, times=None, s=None, readhooks={}):
	''' Compile composites
	
	Parameters
	----------
	qs : list of tuple
	    List of 2-tuples with the entries (vertical level, variable identifier).
	tests : (dict of) list of decider
	    List of deciders defining the composites to be compiled. Composites may be organised 
	    in named groups or may be just a flat list of deciders. Named groups will be saved 
	    together in one file with the different composites as a new dimension.
	times : list of str
	    *Optional*, default ``conf.years``. Times or years to base this composite on.
	s : 2-tuple of int
	    *Optional*, default ``conf.gridsize``. Grid size for the output arrays.
	readhooks : dict of callable
	    *Optional*, default empty. Alternatives used instead of metopen for certain variables.

	Returns
	-------
	dict of np.ndarray with dimensions (bins,y,x)
	    Composite histogram for binnded data for each composite.
	dict of np.ndarray with dimensions (y,x)
	    Composite most frequent value for binned data for each composite.
	dict of np.ndarray with dimensions (y,x)
	    Composite mean for standard data for each composite.
	dict of int
	    Number of contributing time steps for each composite.
	gridlib.grid
	    Metainformation about the data.
	'''
	
	if type(tests) == dict:
		flattened_tests = []
		for grouped_tests in tests.values():
			flattened_tests.extend(grouped_tests)
	else:
		flattened_tests = tests
	
	# Moved here, such that the defaults get evaluated during run time and not during import time
	# At import time, the appropriate context might not have been present yet
	if type(times) == type(None):
		times = conf.years
	if type(s) == type(None):
		s = conf.gridsize
	
	mean = {}
	hist = {}
	mfv  = {}
	cnt  = {}
	for plev, q in qs:
		if q in conf.q_bins:
			for test in flattened_tests:
				hist[test.name,plev,q] = np.zeros((len(conf.q_bins[q]),)+s)
				mfv [test.name,plev,q] = np.zeros(s)
		else:
			for test in flattened_tests:
				mean[test.name,plev,q] = np.zeros(s)
	
	test_qs = set([])
	for test in flattened_tests:
		cnt[test.name] = 0
		if hasattr(test, 'required_qs'):
			test_qs.add(test.required_qs)
	
	test_prev = {}
	test_cur  = {}
	test_next = _get_testdat(times[0], test_qs, readhooks)
	static = None
	for tidx,time in zip(range(len(times)), times):
		test_prev = test_cur
		test_cur  = test_next
		if not time == times[-1]:
			test_next = _get_testdat(times[tidx+1], test_qs, readhooks)
		else: 
			test_next = {}
	
		dat = {}
		plev, q = qs[0]
		if q in readhooks:
			f, dat[plev,q], static_ = readhooks[q](time, plev, q, static=True)
		else:
			f, dat[plev,q], static_ = metopen(conf.file_std % {'time': time, 'plev': plev, 'qf': conf.qf[q]}, q)
		# Inject metadata if missing
		if not hasattr(static_, 't_parsed'):
			static_.t_parsed = metio.str2dts(time)

		for plev, q in qs[1:]:
			if (plev,q) in dat: 
				continue
			if (plev,q) in test_cur:
				dat[plev,q] = test_cur[plev,q]
				continue
			
			if q in readhooks:
				f, dat[plev,q] = readhooks[q](time, plev, q)
			else:
				f, dat[plev,q] = metopen(conf.file_std % {'time': time, 'plev': plev, 'qf': conf.qf[q]}, q, no_static=True)

			if q in LINES:
				dat[plev,q] = utils.mask_lines(dat[plev,q], f[LINES[q]][::])
		
		# Construct static for saving the results, and keep/prepare static_ for the current interval
		if not static:
			static = static_
		else:
			static.t_parsed = np.concatenate((static.t_parsed, static_.t_parsed))
		
		for tidx in range(dat[qs[0]].shape[0]):
			t = static_.t_parsed[tidx]
			for test in flattened_tests:
				if test.match(t, tidx, test_prev, test_cur, test_next):
					for plev, q in qs:
						_add(test.name, plev, q, dat[plev,q][tidx], mean, hist)
					cnt[test.name] += 1

	del dat

	for test in flattened_tests:
		for plev, q in qs:
			if q in conf.q_bins:
				# Save remainder, number of time steps outside the given bins
				hist[test.name,plev,q][-1,::] = cnt[test.name] - hist[test.name,plev,q][:-1,::].sum(axis=0)
				# Calculate the most frequent value based on the histogram
				mfv[test.name,plev,q] = utils.cal_mfv(hist[test.name,plev,q], conf.q_bins[q])
			else:
				mean[test.name,plev,q] /= cnt[test.name]
	
		test.reset()

	return mean, hist, mfv, cnt, static


def save(qs, tests, mean, hist, mfv, cnt, static, s=None):
	''' Save composites using metsave_timeless
	
	Parameters
	----------
	qs : list of tuple
	    List of 2-tuples with the entries (vertical level, variable identifier).
	tests : (dict of) list of decider
	    List of deciders defining the composites to be compiled. Composites may be organised 
	    in named groups or may be just a flat list of deciders. Named groups will be saved 
	    together in one file with the different composites as a new dimension.
	mean : dict of np.ndarray with dimensions (bins,y,x)
	    Composite histogram for binnded data for each composite.
	hist : dict of np.ndarray with dimensions (y,x)
	    Composite most frequent value for binned data for each composite.
	mfv : dict of np.ndarray with dimensions (y,x)
	    Composite mean for standard data for each composite.
	cnt : dict of int
	    Number of contributing time steps for each composite.
	static : gridlib.grid
	    Metainformation about the data.
	s : 2-tuple of int
	    *Optional*, default ``conf.gridsize``. Grid size for the output arrays.
	
	'''

	if type(s) == type(None):
		s = conf.gridsize
	
	if type(tests) == dict:
		for groupname, grouped_tests in tests.items():
			testnames = [test.name for test in grouped_tests]
			
			tosave = {}
			tosave[None,'cnt'] = np.empty((len(grouped_tests),1,1))
			for teidx, test in zip(range(len(grouped_tests)), grouped_tests):
				tosave[None,'cnt'][teidx] = cnt[test.name]

			for plev, q in qs:
				tosave[plev,q] = np.empty((len(grouped_tests),)+s)
				for teidx, test in zip(range(len(grouped_tests)), grouped_tests):
					if q in conf.q_bins:
						tosave[plev,q][teidx,::] = mfv[test.name,plev,q]
						tosave[plev,q+'_hist'][teidx,::] = hist[test.name,plev,q]
					else:
						tosave[plev,q][teidx,::] = mean[test.name,plev,q]

			metsave_timeless(tosave, static, groupname, testnames)
	
	else:
		for test in tests:
			tosave = {}
			for plev, q in qs:
				if q in conf.q_bins:
					tosave[plev,q] = mfv[test.name,plev,q]
					tosave[plev,q+'_hist'] = hist[test.name,plev,q]
				else:
					tosave[plev,q] = mean[test.name,plev,q]

			metsave_timeless(tosave, static, test.name, global_atts={'cnt': cnt[test.name]})
	
	return


# the end
