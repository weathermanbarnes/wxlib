#!/usr/bin/env python 
# -*- encoding: utf-8

raise NotImplementedError('This code is obsolete; '
        'if to be kept within dynlib mid- and long-term, it should '
        'be included as a get_eof function in metio.datasource')

from __future__ import absolute_import, unicode_literals, print_function

from . import utils
from . import tagg
from .shorthands import np, dt, get_instantaneous, get_aggregate, metsave_timeless

import pca_module as pca


# TODO: LINES should be centralised somewhere in the variable definitions!
LINES = {'fronts': 'froff', 'convl': 'cloff', 'defl': 'dloff', 'vorl': 'vloff', 'jetaxis': 'jaoff'}


expvar = ('expvar', None, 'EOF Explained variance', '(0 - 1)')
ecoeff = ('ecoeff', None, 'EOF expansion coefficients', '1')
s.conf.register_variable([expvar, ecoeff, ], [])


def build(qs, eofs, times=None, agg='cal_month', N=10):
	''' Calculate EOFs
	'''

	if type(times) == type(None):
		times = s.conf.years
	
	# TODO: There must be a more general way to do this!
	dates = [dt(min(times),1,1), dt(max(times),12,31,23,59)]
	
	# Fetch and normalise data
	dat = {}
	for plev, q in qs:
		# Try to read aggregated data, otherwise read instantaneous data
		if q in LINES:
			try: 
				datline, static_ = get_aggregate(q, dates, agg, force=True, plevs=plev)
				datloff, none = get_aggregate(LINES[q], dates, agg, force=True, plevs=plev, no_static=True)
				agg_done = True
			except ValueError:
				datline, static_ = get_instantaneous(q, dates, force=True, plevs=plev)
				datloff, none = get_instantaneous(LINES[q], dates, force=True, plevs=plev, no_static=True)
				agg_done = False

			dat_ = utils.smear_lines(datline, datloff)
			del datline, datloff

		else:
			try: 
				dat_, static_ = get_aggregate(q, dates, agg, force=True, plevs=plev)
				agg_done = True
			except ValueError:
				dat_, static_ = get_instantaneous(q, dates, force=True, plevs=plev)
				agg_done = False

		s = dat_.shape
		
		# If only instantaneous data available, aggregate on-the-fly
		if not agg_done:
			aggdates, dat[plev,q] = utils.aggregate(static_.t_parsed, dat_, agg)
			static_ = static_.new_time(aggdates)
		else:
			dat[plev,q] = dat_
		
		# Normalise data
		dat[plev,q] -= dat[plev,q].mean(axis=0)
		dat[plev,q] /= dat[plev,q].std(axis=0)

		# If on lat/lon grid: Take areafactor into account (either as cos- or sqrt(cos)-weighting)
		if static_.gridtype == 'latlon':
			#dat = dat*np.cos(static_.y*np.pi/180.0)
			dat[plev,q] = dat[plev,q] * np.sqrt(np.cos(static_.y*np.pi/180.0))
	
	# Calculate EOFs
	datlin = None
	datlin_mask = {}
	thres = {}
	pattern = {}
	tseries = {}
	static = {}
	for name, tmask, smask in eofs:
		rname = 'rot_' + name

		for plev, q in qs:
			# Make a working copy, apply temporal and spatial masks
			if hasattr(tmask, '__call__'):
				tmask = np.array([tmask(date) for date in static_.t_parsed], dtype=bool)
			dat_ = np.copy(dat[plev,q][tmask,::])
			static[name] = static_.new_time(np.array(static_.t_parsed)[tmask])
			static[rname] = static[name]
			
			if not type(smask) == np.ndarray:
				slc = smask
				smask = np.zeros((s[1], s[2]), dtype=bool)
				smask[slc] = True

			# Flatten to one common EOF input array
			datlin_ = dat_[:,smask]
			datlin_mask[plev,q] = datlin_.std(axis=0) > 0.0
			datlin_ = datlin_[:,datlin_mask[plev,q]]
			if type(datlin) == type(None):
				datlin = datlin_
			else:
				datlin = np.concatenate((datlin, datlin_), axis=1)

		del datlin_, dat_

		print('Starting calculation of EOF `%s`, input shape: %s' % (name, str(datlin.shape)))

		# Do the actual EOF calculations
		tseries[name], pattern_raw, thres[name] = pca.PCA_nipals_c(datlin, PCs=N, standardize=False)
		N_ = len(thres[name])

		# Perform a (normal) varimax rotation
		rpattern_raw, R = utils.varimax(pattern_raw.T)
		tseries[rname] = np.dot(tseries[name], R)
		rtssqsum = (tseries[rname][:,:]**2).sum()
		thres[rname] = np.empty((N_,))
		for i in range(N_):
			# Explained variance is the fraction of the total squared sum of the loadings (= values in the time series)
			thres[rname][i] = (tseries[rname][:,i]**2).sum()/rtssqsum * sum(thres[name])
		
		# unflatten to output data structure
		for plev, q in qs:
			pattern[name,plev,q] = np.zeros((N_, s[1]*s[2]))
			pattern[name,plev,q][:,datlin_mask[plev,q]] = pattern_raw[:,:]
			pattern[name,plev,q] = pattern[name,plev,q].reshape((N_, s[1], s[2]))

			pattern[rname,plev,q] = np.zeros((N_, s[1]*s[2]))
			pattern[rname,plev,q][:,datlin_mask[plev,q]] = rpattern_raw[:,:].T
			pattern[rname,plev,q] = pattern[rname,plev,q].reshape((N_, s[1], s[2]))
	
	return tseries, pattern, thres, static


def save(qs, eofs, tseries, pattern, thres, static):
	''' Save EOFs '''

	vname = ''
	for plev, q in qs:
		vname += '.%s.%s' % (plev, q)

	names = [eof[0] for eof in eofs]
	names.extend(['rot_'+eof[0] for eof in eofs])
	for ename in names:
		name = ename + vname
		if not ename[:4] == 'rot_':
 			idnames = ['EOF %d' % i for i in range(len(thres[ename]))]
		else:
 			idnames = ['REOF %d' % i for i in range(len(thres[ename]))]
		
		tosave = {
			(None, 'expvar'): thres[ename][:,np.newaxis,np.newaxis]
		}
		for plev, q in qs:
			if q in LINES:
				tosave[plev,q+'_freq_pattern'] = pattern[ename,plev,q]
			else:
				tosave[plev,q+'_pattern'] = pattern[ename,plev,q]

		timeseries = {
			(None, 'ecoeff'): tseries[ename],
		}

		metsave_timeless(tosave, static[ename], name, idnames, timeseries=timeseries)

	return



djf_mask = lambda date: date.month in [12,1,2]
mam_mask = lambda date: date.month in [3,4,5]
jja_mask = lambda date: date.month in [6,7,8]
son_mask = lambda date: date.month in [9,10,11]

#
