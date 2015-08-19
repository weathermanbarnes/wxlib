#!/usr/bin/env python
# -*- encoding: utf-8

from ..shorthands import np, metopen
from decider import *



''' This module defines some generally useful deciders to build composites from.

The deciders are based either on fundamental time intervals (day, month, season)
or based on freely available time series for variability indexes. A complete list 
of variability indexes and their abbreviation:

AO
    Arctic Oscillation
NAO 
    North Atlantic Osciallation
EA
    East Atlantic pattern
PNA 
    Pacific-North America pattern
WP
    West Pacific pattern
AAO
    Antarctic Oscillation
MEI
    Multivariate ENSO index
SOI 
    Southern Oscillation Index
MJO
    Madden-Julian Oscillation (index for 8 phases)
'''


def get_months():
	''' Create deciders for each calendar month

	Returns
	-------
	dict of list of decider
	    Deciders for each calendar month, starting with January.
	'''

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

	return {'months' [jan, feb, mar, apr, mai, jun, jul, aug, sep, oct, nov, dec, ], }
	

def get_seasons():
	''' Create deciders for each meteorological season

	Returns
	-------
	dict of list of decider
	    Deciders for each season, starting with DJF.
	'''

	jan, feb, mar, apr, mai, jun, jul, aug, sep, oct, nov, dec = get_months()

	djf = dec | jan | feb
	mam = mar | apr | mai
	jja = jun | jul | aug
	son = sep | oct | nov
	djf.name ='DJF'
	mam.name = 'MAM'
	jja.name = 'JJA'
	son.name = 'SON'

	return {'seasons': [djf, mam, jja, son, ], }

def get_12_seasons():
	''' Create deciders for overlapping 3-month periods centered around each calendar month

	Returns
	-------
	dict of list of decider
	    Deciders for each 3-month period, starting with DJF.
	'''

	jan, feb, mar, apr, mai, jun, jul, aug, sep, oct, nov, dec = get_months()

	djf = dec | jan | feb
	jfm = jan | feb | mar
	fma = feb | mar | apr
	mam = mar | apr | mai
	amj = apr | mai | jun
	mjj = mai | jun | jul
	jja = jun | jul | aug
	jas = jul | aug | sep
	aso = aug | sep | oct
	son = sep | oct | nov
	ond = oct | nov | dec
	ndj = nov | dec | jan

	djf.name = 'DJF'
	jfm.name = 'JFM'
	fma.name = 'FMA'
	mam.name = 'MAM'
	amj.name = 'AMJ'
	mjj.name = 'MJJ'
	jja.name = 'JJA'
	jas.name = 'JAS'
	aso.name = 'ASO'
	son.name = 'SON'
	ond.name = 'OND'
	ndj.name = 'NDJ'

	return {'12seasons': [djf, jfm, fma, mam, amj, mjj, jja, jas, aso, son, ond, ndj], }

def get_diurnal():
	''' Create deciders for daily 6-hour intervals

	Returns
	-------
	dict of list of decider
	    Deciders for each 6-hour interval, starting with [0 UTC, 6 UTC).
	'''

	night   = timeofday('h00-06',  0,  6)
	morning = timeofday('h06-12',  6, 12)
	day     = timeofday('h12-18', 12, 18)
	evening = timeofday('h18-24', 18, 24)

	return {'diurnal_cycle': [night, morning, day, evening], }

def get_nh_daily_indexes():
	''' Create deciders daily northern hemispheric variability indexes

	Currently, there is daily data for NAO and PNA, provided by the `Climate Prediction 
	Center of NOAA <http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/teleconnections.shtml>`_.

	The threshold for the indexes is +/- 1 (standard deviation).

	Returns
	-------
	dict of list of decider
	    Deciders for NAO+, NAO-, PNA+ and PNA-.
	'''

	nao   = metopen('indexes/ts_nao', no_static=True)
	nao_p = ts_lowerbound('NAO+', nao,  1.0)
	nao_n = ts_upperbound('NAO-', nao, -1.0)

	pna   = metopen('indexes/ts_pna', no_static=True)
	pna_p = ts_lowerbound('PNA+', pna,  1.0)
	pna_n = ts_upperbound('PNA-', pna, -1.0)

	return {'NAOd': [nao_p, nao_n], 'PNAd': [pna_p, pna_n, ], }

def get_sh_monthly_indexes():
	''' Create deciders daily northern hemispheric variability indexes

	Currently, there is data for the AAO, provided by the `Climate Prediction 
	Center of NOAA <http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/teleconnections.shtml>`_.

	The threshold for the indexes is +/- 1 (standard deviation).

	Returns
	-------
	dict of list of decider
	    Deciders for AAO+ and AAO-.
	'''

	aao   = metopen('indexes/ts_aao', no_static=True)
	aao_p = ts_lowerbound('AAO+', aao,  1.0)
	aao_n = ts_upperbound('AAO-', aao, -1.0)

	return {'AAO': [aao_p, aao_n, ], }

def get_nh_monthly_indexes():
	''' Create deciders daily northern hemispheric variability indexes

	Currently, there is data for AO, NAO, EA, PNA and WP, all provided by the `Climate Prediction 
	Center of NOAA <http://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml>`_.

	The threshold for the indexes is +/- 1 (standard deviation).

	Returns
	-------
	dict of list of decider
	    Deciders for AO, NAO, EA, PNA and WP, each for positive and negative phases.
	'''

	ao   = metopen('indexes/ts_ao', no_static=True)
	ao_p = ts_lowerbound('AO+', ao,  1.0)
	ao_n = ts_upperbound('AO-', ao, -1.0)
	
	nao   = metopen('indexes/ts_nao_monthly', no_static=True)
	nao_p = ts_lowerbound('NAO+', nao,  1.0)
	nao_n = ts_upperbound('NAO-', nao, -1.0)

	ea    = metopen('indexes/ts_ea', no_static=True)
	ea_p  = ts_lowerbound('EA+', ea,  1.0)
	ea_n  = ts_upperbound('EA-', ea, -1.0)
	
	pna   = metopen('indexes/ts_pna_monthly', no_static=True)
	pna_p = ts_lowerbound('PNA+', pna,  1.0)
	pna_n = ts_upperbound('PNA-', pna, -1.0)

	wp    = metopen('indexes/ts_wp', no_static=True)
	wp_p  = ts_lowerbound('WP+', wp,  1.0)
	wp_n  = ts_upperbound('WP-', wp, -1.0)
	
	return {'AO': [ao_p, ao_n, ], 'NAO': [nao_p, nao_n, ], 'EA': [ea_p, ea_n, ], 
			'PNA': [ pna_p, pna_n, ], 'WP': [wp_p, wp_n, ], }

def get_tropical_indexes():
	''' Create deciders daily northern hemispheric variability indexes

	Currently, there is daily data for MEI and SOI, provided by the `Earth System 
	Research Laboratory <http://www.esrl.noaa.gov/psd/enso/mei/>`_, the Australian 
	`Bureau of Meteorology <http://www.bom.gov.au/climate/current/soi2.shtml>`_ and 
	the `Centre for Australian Weather and Climate Research <http://cawcr.gov.au/staff/mwheeler/maproom/RMM/>`_.

	The threshold for the indexes is +/- 1 (standard deviation).

	Returns
	-------
	dict of list of decider
	    Deciders MEI+, MEI-, SOI+, SOI- and the eight phases of the MJO after Wheeler and Hendon (2004).
	'''

	enso   = metopen('indexes/ts_enso', no_static=True)
	enso_p = ts_lowerbound('MEI+', enso,  1.0)
	enso_n = ts_upperbound('MEI-', enso, -1.0)
	
	soi   = metopen('indexes/ts_soi', no_static=True)
	soi_p = ts_lowerbound('SOI+', soi,  1.0)
	soi_n = ts_upperbound('SOI-', soi, -1.0)

	mjo = metopen('indexes/ts_mjo', no_static=True)
	mjo_1 = ts_equal('MJO_phase1', mjo, 1)
	mjo_2 = ts_equal('MJO_phase2', mjo, 2)
	mjo_3 = ts_equal('MJO_phase3', mjo, 3)
	mjo_4 = ts_equal('MJO_phase4', mjo, 4)
	mjo_5 = ts_equal('MJO_phase5', mjo, 5)
	mjo_6 = ts_equal('MJO_phase6', mjo, 6)
	mjo_7 = ts_equal('MJO_phase7', mjo, 7)
	mjo_8 = ts_equal('MJO_phase8', mjo, 8)

	return {'MEI' : [enso_p, enso_n, ], 'SOI': [soi_p, soi_n, ], 
			'MJO': [mjo_1, mjo_2, mjo_3, mjo_4, mjo_5, mjo_6, mjo_7, mjo_8, ], }

# the end
