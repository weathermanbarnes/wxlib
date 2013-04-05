#!/usr/bin/env python
# -*- encoding: utf-8

from mpl_toolkits.basemap import Basemap


# #############################################################################
# Some often used map projections for separate import
#

# (a) World map
def wmap():
	return Basemap(projection='robin',lon_0=0,resolution='c')
# (b) Northern polar centered map
def npmap():
	return Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
# (c) Southern polar centered map
def spmap():
	return Basemap(projection='spstere',boundinglat=-15,lon_0=0,resolution='l')
# (d) North-Atlantic map
def NAmap():
	return Basemap(projection='lcc', lat_0=55, lat_ts=55, lon_0=-30, resolution='l', 
			width=9000000, height=6000000)
# (e) North-Pacific map
def NPmap():
	return Basemap(projection='lcc', lat_0=50, lat_ts=50, lon_0=-180, resolution='l', 
			width=9000000, height=6000000)
# (f) Australia map
def Ausmap():
	return Basemap(projection='lcc', lat_0=-35, lat_ts=-35, lon_0=120, resolution='l', 
			width=12000000, height=8000000)

# the end
