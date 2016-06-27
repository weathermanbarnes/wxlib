#!/usr/bin/env python
# -*- encoding: utf-8

''' A collection of useful map projections

Much of the research focusses on specific regions. These map projections are intended as 
standard projections covering these regions.
'''

from __future__ import absolute_import, unicode_literals

from mpl_toolkits.basemap import Basemap



# (a) World map
def world():
	''' World map, using the Robin projection

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='robin',lon_0=0,resolution='c', area_thresh=50000)

world.aspect = 1.5

# (b) Northern polar centered map
def n_hemisphere():
	''' Stereographic map, centered on the north pole, covering most of the northern hemisphere

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='npstere',boundinglat=10,lon_0=-50,resolution='c', area_thresh=50000)

n_hemisphere.aspect = 1.0

# (c) Southern polar centered map
def s_hemisphere():
	''' Stereographic map, centered on the south pole, covering most of the southern hemisphere

	Returns
	-------
	Basemap
		map projection instance
	'''
	
	return Basemap(projection='spstere',boundinglat=-10,lon_0=0,resolution='c', area_thresh=50000)

s_hemisphere.aspect = 1.0

# (d) North-Atlantic map
def N_Atlantic():
	''' Map over the North Atlantic, using the Lambert conformal projection

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='lcc', lat_0=55, lat_ts=55, lon_0=-30, resolution='l', 
			width=9000000, height=6000000)

N_Atlantic.aspect = 1.5

# (e) North-Pacific map
def N_Pacific():
	''' Map over the North Pacific, using the Lambert conformal projection

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='lcc', lat_0=50, lat_ts=50, lon_0=-180, resolution='l', 
			width=9000000, height=6000000)

N_Pacific.aspect = 1.5

# (f) Australia map
def Australia():
	''' Map over the Australia, using the Lambert conformal projection

	Returns
	-------
	Basemap
		map projection instance
	'''
	
	return Basemap(projection='lcc', lat_0=-35, lat_ts=-35, lon_0=120, resolution='l', 
			width=12000000, height=8000000)

Australia.aspect = 1.5

# the end
