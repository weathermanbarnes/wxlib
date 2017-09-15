#!/usr/bin/python
# -*- encoding: utf-8

from __future__ import absolute_import, unicode_literals, print_function

from .. import utils, metio
from ..shorthands import np, dt, td, metopen, metsave_timeless, get_static
from ..settings import conf

from scipy.io import savemat
from scipy.interpolate import griddata # interp2d

import mpl_toolkits.basemap
Proj = mpl_toolkits.basemap.pyproj.Proj


# TODO: LINES should be centralised somewhere in the variable definitions!
LINES = {'fronts': 'froff', 'convl': 'cloff', 'defl': 'dloff', 'vorl': 'vloff', 'jetaxis': 'jaoff'}

# TODO: How to avoid hard-coding the rotated grid dimensions here?
PROJGRID_Y, PROJGRID_X = np.meshgrid(np.arange(-1000,1001,40)*1e3, np.arange(-1000,1001,40)*1e3)
PROJGRID_R =  np.sqrt(PROJGRID_X**2 + PROJGRID_Y**2)
S_PROJGRID = PROJGRID_X.shape

cnt = ('cnt', None, 'Number of time steps contributing to the composite', '1')
conf.register_variable([cnt, ], [])


def _add(name, plev, q, dat, mean, hist, x_proj=None, y_proj=None, angle=None):
	''' Add one time step to the composite '''

	if q in conf.q_bins:
		if not type(rot_center) == type(None):
			raise NotImplementedError('Binned variables cannot be rotated yet.')

		for bi in range(len(conf.q_bins[q])-1):
			upper = conf.q_bins[q][bi+1]
			lower = conf.q_bins[q][bi]
			if upper > lower:
				hist[name,plev,q][bi,(dat >= lower).__and__(dat <  upper)] += 1
			else:
				hist[name,plev,q][bi,(dat <  upper).__or__(dat >= lower)] += 1
	else:
		if not type(x_proj) == type(None):
			# Clockwise rotation (following the convention for wind directions)
			x_proj_ =  x_proj*np.cos(angle) + y_proj*np.sin(angle)
			y_proj_ = -x_proj*np.sin(angle) + y_proj*np.cos(angle)
			r_proj_ = np.sqrt(x_proj_**2 + x_proj_**2)

			mask = (r_proj_ <= 1.5e6)

			mean[name,plev,q][:,:] += griddata((x_proj_[mask], y_proj_[mask]), dat[mask], 
					(PROJGRID_Y, PROJGRID_X), method='linear' )
			#interp = interp2d(x_proj_[mask], y_proj_[mask], dat[mask], kind='cubic')
			#mean[name,plev,q][:,:] += interp(PROJGRID_Y[0,:], PROJGRID_X[:,0])

		else:
			mean[name,plev,q][:,:] += dat

	return


def _get_dat(time, test_qs, readhooks, no_static=False):
	''' Retrieve the data used in deciders for testing and for the actual composites '''
	
	static = None
	testdat = {}
	left2request = {}
	for test_plev, test_q in test_qs:
		if not test_q:
			continue
	
		if test_q in readhooks:
			if test_q in left2request:
				left2request[test_q].append(test_plev)
			else:
				left2request[test_q] = [test_plev, ]
		
		elif test_q in LINES:
			testdat[test_plev, test_q] = utils.mask_lines(testdat[test_plev, test_q], f.variables[LINES[test_q]][::])

		else:
			if not no_static and not static:
				f, testdat[test_plev, test_q], static = metopen(conf.file_std % {'time': time, 'plev': test_plev, 'qf': conf.qf[test_q]}, test_q)
			else:
				f, testdat[test_plev, test_q] = metopen(conf.file_std % {'time': time, 'plev': test_plev, 'qf': conf.qf[test_q]}, test_q, no_static=True)

		testdat[test_plev,test_q] = testdat[test_plev,test_q].squeeze()


	# Request all vertical levels of one variable at once for potential more effective read
	for test_q in left2request:
		if not no_static and not static:
			f, dat_, static = readhooks[test_q](time, left2request[test_q], test_q, static=True)
		else:
			f, dat_ = readhooks[test_q](time, left2request[test_q], test_q)

		for pidx, test_plev in zip(range(len(left2request[test_q])), left2request[test_q]):
			testdat[test_plev,test_q] = dat_[:,pidx,:,:].squeeze()

	if no_static:
		return testdat
	else:
		return testdat, static


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

	grid = get_static()
	
	mean = {}
	hist = {}
	mfv  = {}
	cnt  = {}
	for plev, q in qs:
		if q in conf.q_bins:
			for test in flattened_tests:
				if not type(test.rotation_center) == type(None):
					raise NotImplementedError('Binned variables cannot be rotated yet.')
				hist[test.name,plev,q] = np.zeros((len(conf.q_bins[q]),)+s)
				mfv [test.name,plev,q] = np.zeros(s)
		else:
			for test in flattened_tests:
				if not type(test.rotation_center) == type(None):
					mean[test.name,plev,q] = np.zeros(S_PROJGRID)

					# Setup projection and interpolation

					# <ob_tran> is moving the North Pole to a specific given coordinate.
					# However, what we need here is to move a specific location to the 
					# North Pole, which is the inverse transform.
					# They way it's setup here, it should be it's own inverse transform,
					# except for potentially a missing offset in the lon_0.
					rot_proj = Proj(proj='stere', 
						lon_0=test.rotation_center[0], 
						lat_0=test.rotation_center[1], 
					)
					x_proj, y_proj = rot_proj(grid.x, grid.y)

					r_proj = np.sqrt(x_proj**2 + y_proj**2)

					x_proj[np.abs(r_proj) > 1.0e6] = np.nan
					y_proj[np.abs(r_proj) > 1.0e6] = np.nan

					# TODO: Where and how to set angle?
					test.rotargs = (x_proj, y_proj, 0)
				else:
					mean[test.name,plev,q] = np.zeros(s)
					test.rotargs = (None, )*3
	
	test_qs = set([])
	for test in flattened_tests:
		cnt[test.name] = 0
		if hasattr(test, 'required_qs'):
			test_qs.add(test.required_qs)
	
	test_prev = {}
	test_cur  = {}
	test_next = _get_dat(times[0], test_qs, readhooks, no_static=True)
	static = None
	for tidx,time in zip(range(len(times)), times):
		test_prev = test_cur
		test_cur  = test_next
		if not time == times[-1]:
			test_next = _get_dat(times[tidx+1], test_qs, readhooks, no_static=True)
		else: 
			test_next = {}
	
		dat, static_ = _get_dat(time, qs, readhooks)
		# Inject metadata if missing
		if not hasattr(static_, 't_parsed'):
			static_.t_parsed = metio.str2dts(time)

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
						_add(test.name, plev, q, dat[plev,q][tidx], mean, hist, *test.rotargs)
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
			testrots = [type(test.rotation_center) == type(None) for test in grouped_tests]
			if not all(testrots):
				for test in grouped_tests:
					_save_npz(qs, mean, static, test.name, cnt)

			else:
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
			if not type(test.rotation_center) == type(None):
				_save_npz(qs, mean, static, test.name, cnt)

			else:
				tosave = {}
				for plev, q in qs:
					if q in conf.q_bins:
						tosave[plev,q] = mfv[test.name,plev,q]
						tosave[plev,q+'_hist'] = hist[test.name,plev,q]
					else:
						tosave[plev,q] = mean[test.name,plev,q]

				metsave_timeless(tosave, static, test.name, global_atts={'cnt': cnt[test.name]})
	
	return


def _save_npz(qs, mean, static, testname, cnt):
	filename = conf.file_timeless % {'time': metio.dts2str(static.t_parsed), 'name': testname}
	tosave = {}
	for plev, q in qs:
		tosave['%s_%s' % (plev, q)] = mean[testname,plev,q]
	tosave['cnt'] = cnt[testname]
	print('Saving as: %s/%s.npz' % (conf.opath, filename))
	np.savez(conf.opath+'/'+filename+'.npz', **tosave)


# the end
