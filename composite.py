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
plevs = [300, 500, 800] #c.plevs
qs    = ['defabs', 'defang', 'Z', 'u', 'v']

opath = '/work/csp001/composites/minZat300hPa/'

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
		ret.match = lambda date, tidx, field: self.match(date, tidx, field) | b.match(date, tidx, field)
		ret.reset = lambda : self.reset() | b.reset()

		return ret

	def __and__(self, b):
		ret = decider(self.name+'@'+b.name)
		ret.match = lambda date, tidx, field: self.match(date, tidx, field) & b.match(date, tidx, field)
		ret.reset = lambda : self.reset() | b.reset()

		return ret
	
	# The deciding function, returns True/False
	# to be overriden by derived classes.
	def match(self, date, tidx, field):
		return False
	
	# Reset the decider for new vertical level (that is a new pass through the time interval)
	def reset(self):
		return False

class lowerbound_pos(decider):
	# Initialisation
	def __init__(self, name, pos, thres):
		decider.__init__(self, name)
		self.thres = thres
		self.yidx, self.xidx = pos

		return
	
	# Decider
	def match(self, date, tidx, field):
		return field[tidx,self.yidx,self.xidx] >= self.thres
	

class upperbound_pos(lowerbound_pos):
	# Decider
	def match(self, date, tidx, field):
		return field[tidx,self.yidx,self.xidx] < self.thres

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
	def match(self, date, tidx, field):
		while self.dates[self.tidx+1] <= date:
			self.tidx += 1
		
		return self.values[self.tidx] >= self.thres

	# Reset
	def reset(self):
		self.tidx = 0

		return False

class upperbound_ts(lowerbound_ts):
	# Decider
	def match(self, date, tidx, field):
		while self.dates[self.tidx+1] <= date:
			self.tidx += 1
		
		return self.values[self.tidx] < self.thres

class month(decider):
	# Initialisation
	def __init__(self, name, month):
		decider.__init__(self, name)
		self.month = month

		return
	
	# Decider
	def match(self, date, tidx, field):
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
djf.name ='DJF'
mam.name = 'MAM'
jja.name = 'JJA'
son.name = 'SON'

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

alps     = upperbound_pos('Alps_TB',      (88, 371), 89385.6010909)
bagheran = upperbound_pos('Bagheran_TB', (114, 480), 93406.560068)
greenlnd = upperbound_pos('Greenland_TB', (51, 278), 86043.961205)
rockies  = upperbound_pos('Rockies_TB',  (114, 480), 93406.560068)

pacsec30 = upperbound_pos('PacSec_30N',  (120,   0), 93321.2289874)
pacsec35 = upperbound_pos('PacSec_35N',  (110,   0), 92138.791619)
pacsec40 = upperbound_pos('PacSec_40N',  (100,   0), 89810.2522804)
pacsec45 = upperbound_pos('PacSec_45N',   (90,   0), 87825.1021409)
pacsec50 = upperbound_pos('PacSec_50N',   (80,   0), 86791.546998)
pacsec55 = upperbound_pos('PacSec_55N',   (70,   0), 86541.7791012)

atlsec35 = upperbound_pos('AtlSec_35N',  (110, 300), 92042.3654283)
atlsec40 = upperbound_pos('AtlSec_40N',  (100, 300), 90781.3529516)
atlsec45 = upperbound_pos('AtlSec_45N',   (90, 300), 89316.0434635)
atlsec50 = upperbound_pos('AtlSec_50N',   (80, 300), 87986.0890262)
atlsec55 = upperbound_pos('AtlSec_55N',   (70, 300), 87042.1297982)
atlsec60 = upperbound_pos('AtlSec_60N',   (60, 300), 86423.7978621)

sibsec45 = upperbound_pos('SibSec_45N',   (90, 510), 90659.7261204)
sibsec50 = upperbound_pos('SibSec_50N',   (80, 510), 89170.244846)
sibsec55 = upperbound_pos('SibSec_55N',   (70, 510), 87829.8385007)
sibsec60 = upperbound_pos('SibSec_60N',   (60, 510), 86805.0279739)
sibsec65 = upperbound_pos('SibSec_65N',   (50, 510), 85983.7165583)
sibsec70 = upperbound_pos('SibSec_70N',   (40, 510), 85582.3673396)

aussec35 = upperbound_pos('AusSec_35S',  (250, 600), 90282.5965711)
aussec40 = upperbound_pos('AusSec_40S',  (260, 600), 87925.7773702)
aussec45 = upperbound_pos('AusSec_45S',  (270, 600), 85857.7848112)
aussec50 = upperbound_pos('AusSec_50S',  (280, 600), 84169.2246727)
aussec55 = upperbound_pos('AusSec_55S',  (290, 600), 83053.8476946)
aussec60 = upperbound_pos('AusSec_60S',  (300, 600), 82323.6388425)
aussec65 = upperbound_pos('AusSec_65S',  (310, 600), 81855.9116357)

#tests = [jan, feb, mar, apr, mai, jun, jul, aug, sep, oct, nov, dec, ]
#tests = [djf, mam, jja, son, ]
#tests = [ao_p & djf, ao_n & djf, nao_p & djf, nao_n & djf, aao_p & djf, aao_n & djf, 
#	 pna_p & djf, pna_n & djf, enso_p & djf, enso_n & djf, ]
#tests = [ao_p & jja, ao_n & jja, nao_p & jja, nao_n & jja, aao_p & jja, aao_n & jja, 
#	 pna_p & jja, pna_n & jja, enso_p & jja, enso_n & jja, ]
#tests = [alps & djf, bagheran & djf, greenlnd & djf, rockies & djf,
#	 atlsec35 & djf, atlsec40 & djf, atlsec45 & djf, atlsec50 & djf, atlsec55 & djf, atlsec60 & djf, 
#	 sibsec45 & djf, sibsec50 & djf, sibsec55 & djf, sibsec60 & djf, sibsec65 & djf, sibsec70 & djf,
#	 aussec35 & jja, aussec40 & jja, aussec45 & jja, aussec50 & jja, aussec55 & jja, aussec60 & jja, aussec65 & jja,
#	 pacsec30 & djf, pacsec35 & djf, pacsec40 & djf, pacsec45 & djf, pacsec50 & djf, pacsec55 & djf,  ]
tests = [alps & jja, bagheran & jja, greenlnd & jja, rockies & jja, 
	 atlsec35 & jja, atlsec40 & jja, atlsec45 & jja, atlsec50 & jja, atlsec55 & jja, atlsec60 & jja,
 	 sibsec45 & jja, sibsec50 & jja, sibsec55 & jja, sibsec60 & jja, sibsec65 & jja, sibsec70 & jja,
	 aussec35 & djf, aussec40 & djf, aussec45 & djf, aussec50 & djf, aussec55 & djf, aussec60 & djf, aussec65 & djf,
	 pacsec30 & jja, pacsec35 & jja, pacsec40 & jja, pacsec45 & jja, pacsec50 & jja, pacsec55 & jja,  ]


test_q    = 'Z'
test_plev = 300

# ---------------------------------------------------------------------
# Building the composites

print 'Preparing'

f, oro = metopen('static', 'oro', cut=c.std_slice[1:])
s = oro.shape
f.close()
del oro

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
			mfv[j,i] = (bins[bi+1]+bins[bi])/2.0

	return mfv

for plev in plevs:
	dat = {}
	mean = {}
	hist = {}
	mfv  = {}
	cnt  = np.zeros((len(tests),))
	for q in qs:
		if q in c.bins:
			hist[q] = np.zeros((len(tests), len(c.bins[q]), s[0], s[1]))
			mfv [q] = np.zeros((len(tests), s[0], s[1]))
		else:
			mean[q] = np.zeros((len(tests), s[0], s[1]))

	for yr in years:
		f, testdat = metopen(c.file_std % (yr, test_plev, test_q), c.q[test_q])
	
		print plev, yr

		for q in qs:
			f, dat[q] = metopen(c.file_std % (yr, plev, q), c.q[q])
		
		t0 = dt(yr, 1, 1)
		for tidx in range(testdat.shape[0]):
			t = t0 + tidx*td(0.25,0)
			for ti in range(len(tests)):
				if tests[ti].match(t, tidx, testdat):
					for q in qs:
						add(ti, q, dat[q][tidx])
					cnt[ti] += 1

	print plev, 'Postprocessing', cnt

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

		tosave['cnt'] = cnt[ti]
		
		np.savez(opath+'%s_composite.%s.npz' % (tests[ti].name, plev), **tosave)
		
		tests[ti].reset()

# the end
