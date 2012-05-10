#!/usr/bin/python
# -*- encoding: utf-8

import sys
import numpy as np
import static as c
from metopen import metopen
from datetime import datetime as dt, timedelta as td

# ---------------------------------------------------------------------
# Data base

years = c.years
plevs = c.plevs
qs    = ['defabs', 'defang', 'Z']

# ---------------------------------------------------------------------
# Defining possible testers

class decider(object):
	# Constructor; potentially overriden by derived classes
	def __init__(self, name):
		self.name = name

		return
	
	# Construct new deciders by logically connecting two existing decider instances
	def __or__(self, b):
		ret = decider(self.name+b.name)
		ret.match = lambda date, field: self.match(date, field) | b.match(date, field)

		return ret

	def __and__(self, b):
		ret = decider(self.name+'@'+b.name)
		ret.match = lambda date, field: self.match(date, field) & b.match(date, field)

		return ret
	
	# The deciding function, returns True/False
	# to be overriden by derived classes.
	def match(self, date, field):
		return False

class lowerbound_pos(decider):
	# Initialisation
	def __init__(self, name, pos, thres):
		decider.__init__(self, name)
		self.thres = thres
		self.yidx, self.xidx = pos

		return
	
	# Decider
	def match(self, date, field):
		return field[::,self.yidx,self.xidx] >= self.thres

class upperbound_pos(lowerbound_pos):
	# Decider
	def match(self, date, field):
		return field[::,self.yidx,self.xidx] < self.thres

class lowerbound_ts(decider):
	# Initialisation
	def __init__(self, name, ts, thres):
		decider.__init__(self, name)
		self.tidx   = 0
		self.dates  = ts['dates']
		self.values = ts['values']
		self.thres  = thres

		return
	
	# Decider
	def match(self, date, field):
		while self.dates[self.tidx] < date:
			self.tidx += 1
		
		return self.values[self.tidx] >= self.thres

class upperbound_ts(lowerbound_ts):
	# Decider
	def match(self, date, field):
		while self.dates[self.tidx] < date:
			self.tidx += 1
		
		return self.values[self.tidx] < self.thres

class month(decider):
	# Initialisation
	def __init__(self, name, month):
		decider.__init__(self, name)
		self.month = month

		return
	
	# Decider
	def match(self, date, field):
		return date.month == self.month

# ---------------------------------------------------------------------
# Building the testers

jan = month('Jan', 1)
feb = month('Feb', 2)
mar = month('Mar', 3)
apr = month('Apr', 4)
mai = month('Mai', 5)
jun = month('Jun', 6)
jul = month('Jul', 7)
aug = month('Aug', 8)
sep = month('Sep', 9)
oct = month('Oct', 10)
nov = month('Nov', 11)
dec = month('Dec', 12)

djf = dec | jan | feb
mam = mar | apr | mai
jja = jun | jul | aug
son = sep | oct | nov

ao   = np.load('ts_ao.npz')
ao_p = lowerbound_ts('AO+', ao,  1.0)
ao_n = upperbound_ts('AO-', ao, -1.0)

nao   = np.load('ts_nao.npz')
nao_p = lowerbound_ts('NAO+', nao,  1.0)
nao_n = upperbound_ts('NAO-', nao, -1.0)

aao   = np.load('ts_aao.npz')
aao_p = lowerbound_ts('AAO+', aao,  1.0)
aao_n = upperbound_ts('AAO-', aao, -1.0)

pna   = np.load('ts_pna.npz')
pna_p = lowerbound_ts('PNA+', pna,  1.0)
pna_n = upperbound_ts('PNA-', pna, -1.0)

enso   = np.load('ts_enso.npz')
enso_p = lowerbound_ts('ENSO+', enso,  1.0)
enso_n = upperbound_ts('ENSO-', enso, -1.0)

#tests = [jan, feb, mar, apr, mai, jun, jul, aug, sep, oct, nov, dec, djf, mam, jja, son ]
tests = [ao_p, ao_n, nao_p, nao_n, aao_p, aao_n, pna_p, pna_n, enso_p, enso_n]

test_q    = 'defabs'
test_plev = 700

# ---------------------------------------------------------------------
# Building the composites

print 'Preparing'

f, oro = metopen('static', 'oro', cut=c.std_slice[1:])
s = oro.shape
f.close()
del oro

mean = {}
hist = {}
mfv  = {}
cnt  = np.zeros((len(tests),))
for q in qs:
	if q in c.bins:
		hist[q] = np.zeros((len(tests), len(c.bins[q]), s[0], s[1]))
		mfv [q] = np.zeros((len(tests), s[0], s[1]), dtype='i4')
	else:
		mean[q] = np.zeros((len(tests), s[0], s[1]))

def add(ti, q, dat):
	if q in c.bins:
		for bi in range(len(c.bins[q])-1):
			upper = c.bins[q][bi+1]
			lower = c.bins[q][bi]
			if upper > lower:
				hist[q][ti,bi,(dat >= lower).__and__(dat <  upper)] += 1
			else:
				hist[q][ti,bi,(dat <  upper).__or__(dat >= lower)] += 1
	else:
		mean[q][ti,:,:] += dat

	return

def cal_mfv(hist, bins):
	mfv = np.zeros(s)
	for j in range(s[0]):
		for i in range(s[1]):
			bi = hist[:,j,i].argmax()
			mfv = (bins[bi+1]+bins[bi])/2.0

	return mfv

for plev in plevs:
	dat = {}

	for yr in years:
		f, testdat = metopen(c.file_std % (yr, test_plev, test_q), c.q[test_q])
	
		print plev, yr

		for q in qs:
			f, dat[q] = metopen(c.file_std % (yr, plev, q), c.q[q])
		
		t0 = dt(yr, 1, 1)
		for tidx in range(testdat.shape[0]):
			t = t0 + tidx*td(0.25,0)
			for ti in range(len(tests)):
				if tests[ti].match(t, testdat):
					for q in qs:
						add(ti, q, dat[q][tidx])
					cnt[ti] += 1

	print plev, 'Postprocessing'

	for ti in range(len(tests)):
		for q in qs:
			if q in c.bins:
				hist[q][ti,-1] = cnt[ti] - hist[q][ti,:-1].sum(axis=0)
				mfv [q][ti]    = cal_mfv(hist[q][ti], c.bins[q])
			else:
				mean[q][ti] /= cnt[ti]

	# ---------------------------------------------------------------------
	# Saving the composites

	for ti in range(len(tests)):
		tosave = {}
		for q in qs:
			if q in c.bins:
				tosave['%s_hist' % q] = hist[q][ti]
				tosave['%s_mfv' % q]  =  mfv[q][ti]
			else:
				tosave['%s_mean' % q] = mean[q][ti]
		
		np.savez('%s_composite.%d.npz' % (tests[ti].name, plev), **tosave)

# the end
