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
qs    = ['defabs', 'defang', 'Z', 'u', 'v']

opath = '/work/csp001/composites/'

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

		return ret

	def __and__(self, b):
		ret = decider(self.name+'@'+b.name)
		ret.match = lambda date, tidx, field: self.match(date, tidx, field) & b.match(date, tidx, field)

		return ret
	
	# The deciding function, returns True/False
	# to be overriden by derived classes.
	def match(self, date, tidx, field):
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
		while self.dates[self.tidx] < date:
			self.tidx += 1
		
		return self.values[self.tidx] >= self.thres

class upperbound_ts(lowerbound_ts):
	# Decider
	def match(self, date, tidx, field):
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

alps     = lowerbound_pos('Alps_TB',      (88, 371), 0.000237193031353)
bagheran = lowerbound_pos('Bagheran_TB', (114, 480), 0.000134757516207)
greenlnd = lowerbound_pos('Greenland_TB', (51, 278), 0.000224143630476)
rockies  = lowerbound_pos('Rockies_TB',  (114, 480), 0.000134757516207)

pacsec30 = lowerbound_pos('PacSec_30N',  (120,   0), 0.00014751413255)
pacsec35 = lowerbound_pos('PacSec_35N',  (110,   0), 0.000172332831426)
pacsec40 = lowerbound_pos('PacSec_40N',  (100,   0), 0.000191705010366)
pacsec45 = lowerbound_pos('PacSec_45N',   (90,   0), 0.000210880651139)
pacsec50 = lowerbound_pos('PacSec_50N',   (80,   0), 0.000217700013309)
pacsec55 = lowerbound_pos('PacSec_55N',   (70,   0), 0.000212105835089)

atlsec35 = lowerbound_pos('AtlSec_35N',  (110, 300), 0.000181586903636)
atlsec40 = lowerbound_pos('AtlSec_40N',  (100, 300), 0.00020357221365)
atlsec45 = lowerbound_pos('AtlSec_45N',   (90, 300), 0.000218390370719)
atlsec50 = lowerbound_pos('AtlSec_50N',   (80, 300), 0.000240608642343)
atlsec55 = lowerbound_pos('AtlSec_55N',   (70, 300), 0.000243240370764)
atlsec60 = lowerbound_pos('AtlSec_60N',   (60, 300), 0.000235398852965)

sibsec45 = lowerbound_pos('SibSec_45N',   (90, 510), 0.00019151752349)
sibsec50 = lowerbound_pos('SibSec_50N',   (80, 510), 0.000185327575309)
sibsec55 = lowerbound_pos('SibSec_55N',   (70, 510), 0.000183245167136)
sibsec60 = lowerbound_pos('SibSec_60N',   (60, 510), 0.000193859974388)
sibsec65 = lowerbound_pos('SibSec_65N',   (50, 510), 0.000188957434148)
sibsec70 = lowerbound_pos('SibSec_70N',   (40, 510), 0.000188765014173)

aussec35 = lowerbound_pos('AusSec_35S',  (250, 600), 0.000177936613909)
aussec40 = lowerbound_pos('AusSec_40S',  (260, 600), 0.000183534648386)
aussec45 = lowerbound_pos('AusSec_45S',  (270, 600), 0.000195973771042)
aussec50 = lowerbound_pos('AusSec_50S',  (280, 600), 0.000204542753636)
aussec55 = lowerbound_pos('AusSec_55S',  (290, 600), 0.000203898554901)
aussec60 = lowerbound_pos('AusSec_60S',  (300, 600), 0.000202376424568)
aussec65 = lowerbound_pos('AusSec_65S',  (310, 600), 0.000192796826013)

#tests = [jan, feb, mar, apr, mai, jun, jul, aug, sep, oct, nov, dec, ]
#tests = [djf, mam, jja, son, ]
#tests = [ao_p & djf, ao_n & djf, nao_p & djf, nao_n & djf, aao_p & djf, aao_n & djf, 
#	 pna_p & djf, pna_n & djf, enso_p & djf, enso_n & djf, ]
#tests = [ao_p & jja, ao_n & jja, nao_p & jja, nao_n & jja, aao_p & jja, aao_n & jja, 
#	 pna_p & jja, pna_n & jja, enso_p & jja, enso_n & jja, ]
#tests = [alps & djf, bagheran & djf, greenlnd & djf, rockies & djf,
#	 alps & jja, bagheran & jja, greenlnd & jja, rockies & jja, ]
#tests = [pacsec30 & djf, pacsec35 & djf, pacsec40 & djf, pacsec45 & djf, pacsec50 & djf, pacsec55 & djf,
#	 pacsec30 & jja, pacsec35 & jja, pacsec40 & jja, pacsec45 & jja, pacsec50 & jja, pacsec55 & jja, ]
#tests = [atlsec35 & djf, atlsec40 & djf, atlsec45 & djf, atlsec50 & djf, atlsec55 & djf, atlsec60 & djf, 
#	 atlsec35 & jja, atlsec40 & jja, atlsec45 & jja, atlsec50 & jja, atlsec55 & jja, atlsec60 & jja, ]
#tests = [sibsec45 & djf, sibsec50 & djf, sibsec55 & djf, sibsec60 & djf, sibsec65 & djf, sibsec70 & djf,
#	 sibsec45 & jja, sibsec50 & jja, sibsec55 & jja, sibsec60 & jja, sibsec65 & jja, sibsec70 & jja, ]
#tests = [aussec35 & djf, aussec40 & djf, aussec45 & djf, aussec50 & djf, aussec55 & djf, aussec60 & djf, aussec65 & djf,
#	 aussec35 & jja, aussec40 & jja, aussec45 & jja, aussec50 & jja, aussec55 & jja, aussec60 & jja, aussec65 & jja ]

test_q    = 'defabs'
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
			mfv = (bins[bi+1]+bins[bi])/2.0

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
			mfv [q] = np.zeros((len(tests), s[0], s[1]), dtype='i4')
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

# the end
