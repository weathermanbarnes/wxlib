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

world.aspect = 2.0

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

# (c) Northern polar centered map, focussing on extratropics
def n_extratropics():
	''' Stereographic map, centered on the north pole, covering most of the northern hemisphere

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='npstere',boundinglat=35,lon_0=-50,resolution='c', area_thresh=10000)

n_extratropics.aspect = 1.0

# (d) Northern polar centered map, focussing on Nordic and Polar Seas
def atlantic_arctic():
	''' Stereographic map, centered on the north pole, covering most of the northern hemisphere

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='npstere',boundinglat=53,lon_0=0,resolution='l', area_thresh=10000)

atlantic_arctic.aspect = 1.0

# (e) Southern polar centered map
def s_hemisphere():
	''' Stereographic map, centered on the south pole, covering most of the southern hemisphere

	Returns
	-------
	Basemap
		map projection instance
	'''
	
	return Basemap(projection='spstere',boundinglat=-10,lon_0=0,resolution='c', area_thresh=50000)

s_hemisphere.aspect = 1.0

# (f) Southern polar centered map, focussing on extratropics
def s_extratropics():
	''' Stereographic map, centered on the north pole, covering most of the northern hemisphere

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='spstere',boundinglat=-35,lon_0=0,resolution='c', area_thresh=10000)

s_extratropics.aspect = 1.0

# (g) North-Atlantic map
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

# (h) Polar North Atlantic: Fram-strait to Greenland-Iceland-Scotland ridge
def Polar_N_Atlantic():
	''' Map over the North Atlantic, using the Lambert conformal projection

	Returns
	-------
	Basemap
		map projection instance
	'''

	return Basemap(projection='lcc', lat_0=67, lat_ts=67, lon_0=-30, resolution='l', 
			width=5000000, height=5000000)

Polar_N_Atlantic.aspect = 1.0

# (i) North-Pacific map
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

# (j) Australia map
def Australia():
	''' Map over Australia, using the Lambert conformal projection

	Returns
	-------
	Basemap
		map projection instance
	'''
	
	return Basemap(projection='lcc', lat_0=-35, lat_ts=-35, lon_0=120, resolution='l', 
			width=12000000, height=8000000)

Australia.aspect = 1.5

# (g) Maps from iveret
# -> Immediate surroundings of Bergen
def Bergen():
	return Basemap(projection='lcc', lat_0=60, lat_ts=60, lon_0=15, resolution='h', 
			llcrnrlon=4.75, llcrnrlat=60.1, urcrnrlon=5.8, urcrnrlat=60.8)
Bergen.aspect = 1.0

# -> Approximately Vestland fylke
def Vestland():
	return Basemap(projection='lcc', lat_0=60, lat_ts=60, lon_0=15, resolution='h', 
			llcrnrlon=2.0, llcrnrlat=58.5, urcrnrlon=7.97, urcrnrlat=62.5)
Vestland.aspect = 1.0

# -> Southern Norway
def S_Norway():
	return Basemap(projection='lcc', lat_0=60, lat_ts=60, lon_0=15, resolution='i', 
			llcrnrlon=-1.5, llcrnrlat=55.4, urcrnrlon=13, urcrnrlat=65.0)
S_Norway = 1.0

# -> Scandinavia, covering about the entire AROME MetCoop domain
def Scandinavia():
	return Basemap(projection='lcc', lat_0=63.3, lat_ts=63.3, lon_0=15, resolution='i', 
			llcrnrlon=0.3, llcrnrlat=50.3, urcrnrlon=54.2, urcrnrlat=71.5)
Scandinavia.aspect = 1.125

# -> North Atlantic
def N_Atlantic_forecast():
	return Basemap(projection='lcc', lat_0=58, lat_ts=58, lon_0=-30, resolution='l', 
			width=10800000, height=7200000)
N_Atlantic_forecast.aspect = 1.5

# -> Barents Sea, covering about the entire AROME Arctic domain
def Barents_Sea():
	return Basemap(projection='lcc', lat_0=77.5, lat_ts=77.5, lon_0=-25, resolution='i', 
			llcrnrlon=-16.0, llcrnrlat=69.5, urcrnrlon=67.5, urcrnrlat=71.5)
Barents_Sea.aspect = 1.25


# the end
