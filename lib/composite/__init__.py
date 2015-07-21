#!/usr/bin/python
# -*- encoding: utf-8

import dynlib.utils 
from dynlib.shorthands import np, dt, td, metopen
from dynlib.settings import conf

from scipy.io import savemat


def _add(name, plev, q, dat, mean, hist):
	''' Add one time step to the composite '''

	if q in conf.bins:
		for bi in range(len(conf.bins[q])-1):
			upper = conf.bins[q][bi+1]
			lower = conf.bins[q][bi]
			if upper > lower:
				hist[name,plev,q][bi,(dat >= lower).__and__(dat <  upper)] += 1
			else:
				hist[name,plev,q][bi,(dat <  upper).__or__(dat >= lower)] += 1
	else:
		mean[name,plev,q][:,:] += dat

	return


def _get_testdat(yr, test_qs):
	''' Retrieve the data used in deciders for testing '''

	testdat = {}
	for test_plev, test_q in test_qs:
		if not test_q:
			continue
	
		if test_q[-6:] == 'front':
			if (test_plev, 'front') in testdat:
				continue
			f, fronts = metopen(conf.file_std % {'time': yr, 'plev': test_plev, 'q': conf.qi['front']}, 'front', no_static=True)
			c, w, s = dynlib.utils.mask_fronts(fronts, f['froff'])
			testdat[test_plev, 'cold_front', ] = c
			testdat[test_plev, 'warm_front'] = w
			testdat[test_plev, 'stat_front'] = s
			testdat[test_plev, 'front'] = np.logical_or(c, np.logical_or(w, s))
		elif test_q == 'convl':
			f, convls = metopen(conf.file_std % {'time': yr, 'plev': test_plev, 'q': conf.qi['convl']}, 'convl', no_static=True)
			c = dynlib.utils.mask_lines(convls, f['cloff'])
			testdat[test_plev, 'convl'] = c
		elif test_q == 'defl':
			f, defls = metopen(conf.file_std % {'time': yr, 'plev': test_plev, 'q': conf.qi['defl']}, 'defl', no_static=True)
			d = dynlib.utils.mask_lines(defls, f['dloff'])
			testdat[test_plev, 'defl'] = d
		elif test_q == 'vorl':
			f, vorls = metopen(conf.file_std % {'time': yr, 'plev': test_plev, 'q': conf.qi['vorl']}, 'vorl', no_static=True)
			v = dynlib.utils.mask_lines(vorls, f['vloff'])
			testdat[test_plev, 'vorl'] = v
		elif test_q == 'jetaxis':
			f, jaxis = metopen(conf.file_std % {'time': yr, 'plev': test_plev, 'q': conf.qi['jetaxis']}, 'jetaxis', no_static=True)
			v = dynlib.utils.mask_lines(f.variables['jetaxis'], f.variables['jaoff'])
			testdat[test_plev, 'jetaxis'] = v
		else:
			f, testdat[test_plev, test_q] = metopen(conf.file_std % {'time': yr, 'plev': test_plev, 'q': conf.qi[test_q]}, test_q, no_static=True)
	
	return testdat


def build(qs, tests, years=conf.years, s=conf.gridsize):
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
	'''
	
	# TODO: LINES should be centralised somewhere in the variable definitions!
	LINES = {'fronts': 'froff', 'convl': 'cloff', 'defl': 'dloff', 'vorl': 'vloff', 'jetaxis': 'jaoff'}

	if type(tests) == dict:
		flattened_tests = []
		for grouped_tests in tests.values():
			flattened_tests.extend(grouped_tests)
	else:
		flattened_tests = tests
	
	mean = {}
	hist = {}
	mfv  = {}
	cnt  = {}
	for plev, q in qs:
		if (plev, q) in conf.bins:
			for test in flattened_tests:
				hist[test.name,plev,q] = np.zeros((len(conf.bins[q]),)+s)
				mfv [test.name,plev,q] = np.zeros(s)
		else:
			for test in flattened_tests:
				mean[test.name,plev,q] = np.zeros(s)
	
	test_qs = []
	for test in flattened_tests:
		test_qs.extend(test.required_qs)
	
	test_prev = {}
	test_cur  = {}
	test_next = _get_testdat(years[0], test_qs)
	for yr in years:
		test_prev = test_cur
		test_cur  = test_next
		if not yr == years[-1]:
			test_next = _get_testdat(yr+1, test_qs)
		else: 
			test_next = {}
	
		dat = {}
		for plev, q in qs:
			if (plev,q) in dat: continue
			
			f, dat[plev,q] = metopen(conf.file_std % {'time': yr, 'plev': plev, 'qf': conf.qf[q]}, q, no_static=True)

			if q in LINES:
				dat[plev,q] = dynlib.utils.mask_lines(dat[plev,q], f[LINES[q]][::])
		
		t0 = dt(yr, 1, 1)
		for tidx in range(dat[qs[0]].shape[0]):
			t = t0 + tidx*td(0.25,0)
			for test in flattened_tests:
				if test.match(t, tidx, test_prev, test_cur, test_next):
					for plev, q in qs:
						_add(test.name, plev, q, dat[plev,q][tidx], mean, hist)
					cnt[test.name] += 1

	del dat

	for test in flattened_tests:
		for plev, q in qs:
			if q in conf.bins:
				# Save remainder, number of time steps outside the given bins
				hist[test.name,plev,q][-1,::] = cnt[test.name] - hist[test.name,plev,q][:-1,::].sum(axis=0)
				# Calculate the most frequent value based on the histogram
				mfv[test.name,plev,q] = dynlib.utils.cal_mfv(hist[test.name,plev,q], conf.bins[q])
			else:
				mean[test.name,plev,q] /= cnt[test.name]
	
		test.reset()

	return mean, hist, mfv, cnt


def save(qs, tests, mean, hist, mfv, cnt):
	# TODO: Where to get static from?
	if type(tests) == dict:
		for groupname, grouped_tests in tests:
			testnames = [test.name for test in grouped_tests]
			
			tosave = {}
			tosave[None,'cnt'] = np.empty((len(grouped_tests),1,1))
			for teidx, test in zip(range(len(grouped_tests)), grouped_tests):
				tosave[None,'cnt'][teidx] = cnt[test.name]

			for plev, q in qs:
				tosave[plev,q] = np.empty((len(grouped_tests),)+s)
				for teidx, test in zip(range(len(grouped_tests)), grouped_tests):
					if q in conf.bins:
						tosave[plev,q][teidx,::] = mfv[test.name,plev,q]
						tosave[plev,q+'_hist'][teidx,::] = hist[test.name,plev,q]
					else:
						tosave[plev,q][teidx,::] = mean[test.name,plev,q]

			metsave_timeless(tosave, static, groupname, testnames)
	
	else:
		for test in tests:
			tosave = {}
			for plev, q in qs:
				if q in conf.bins:
					tosave[plev,q] = mfv[test.name,plev,q]
					tosave[plev,q+'_hist'] = hist[test.name,plev,q]
				else:
					tosave[plev,q] = mean[test.name,plev,q]

			metsave_timeless(tosave, static, test.name, global_atts={'cnt': cnt[test.name]})
	
	return


# the end
