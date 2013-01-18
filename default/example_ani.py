#!/usr/bin/python
# -*- encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, cm as bmcm
from metopen import metopen
from utils import concat1

from datetime import timedelta as td

import static as c


import figures as f

lat = f.lat
lon = f.lon
oro = f.oro
orolevs = f.orolevs
deflevs = np.arange(4.5,27.1,1.5)

# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def map_Q(date_start, date_end, q='defabs', plev=700, cmap=None):
	date_int = [date_start, date_end]
	dat = f._get_instantaneous(q, date_int, plevs=plev, tavg=False)
	daZ = f._get_instantaneous('Z', date_int, plevs=plev)

	dat[:,daZ < oro] = np.nan

	fig = plt.figure()

	def update(t):
		plt.clf()

		m = Basemap(projection='npstere',boundinglat=30,lon_0=-50,resolution='l')
		x,y = m(lon,lat)
		m.drawcoastlines()
		m.contour(x,y, oro, range(5000,80001,5000), colors='k', alpha=0.4, zorder=2)
		m.contourf(x, y, dat[t,:,:], 25, cmap=cmap, zorder=1)
		m.drawparallels(range(15,80,5))
		m.drawmeridians(range(0,360,30))
		plt.title(str(date_start+t*td(0.25)))
		plt.colorbar()
	
	ani = animation.FuncAnimation(fig, update, interval=500, frames=range(dat.shape[0]))
	plt.show()

	return


# Contour map of averaged wind on top of a oro map.
def map_barb(date_start, date_end, plev=700, quiver=False):
	date_int = [date_start, date_end]
	defabs = f._get_instantaneous('defabs', date_int, plevs=plev, tavg=False)
	dau = f._get_instantaneous('u', date_int, plevs=plev, tavg=False)
	dav = f._get_instantaneous('v', date_int, plevs=plev, tavg=False)
	daZ = f._get_instantaneous('Z', date_int, plevs=plev)

	dau[:,daZ < oro] = np.nan
	dav[:,daZ < oro] = np.nan
	
	fig = plt.figure(figsize=(14.4,10.8), dpi=72)
	
	def update(t):
		plt.clf()

		m = Basemap(projection='npstere',boundinglat=30,lon_0=-50,resolution='l')
		x,y = m(lon,lat)
		ut,vt,xt,yt = m.transform_vector(dau[t,::-1,:],dav[t,::-1,:],lon[0,:],lat[::-1,0],80,80,returnxy=True)
		m.drawcoastlines()
		m.contour(x,y, oro, range(5000,80001,5000), colors='k', alpha=0.4, zorder=2)
		pabs = m.contourf(x, y, defabs[t,:,:]*86400, np.arange(31), zorder=1)
		if not quiver:
			m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5, zorder=3)
		else:
			m.quiver(xt, yt, ut, vt, zorder=3)
		m.drawparallels(range(15,80,5))
		m.drawmeridians(range(0,360,30))
		plt.title(str(date_start+t*td(0.25)))
		plt.colorbar()
	
	ani = animation.FuncAnimation(fig, update, interval=500, frames=range(dau.shape[0]))
	ani.save('test.mp4', clear_temp=False, frame_prefix='_tmq_db')
	plt.show()

	return


# Contour map of averaged deformation vector on top of a oro map.
def map_deform(date_start, date_end, plev=700, img_prefix='_tmp_da'):
	date_int = [date_start, date_end]
	defabs = f._get_instantaneous('defabs', date_int, plevs=plev, tavg=False)
	defang = f._get_instantaneous('defang', date_int, plevs=plev, tavg=False)
	daZ    = f._get_instantaneous('Z', date_int, plevs=plev)

	daZ = daZ[::-1,:]
	print defabs.shape, defang.shape, daZ.shape
	
	defdex = np.cos(defang[:,:]) *defabs
	defdey = np.sin(defang[:,:]) *defabs
	
	defabs[:,daZ < oro] = np.nan
	defdex[:,daZ < oro] = np.nan
	defdey[:,daZ < oro] = np.nan
	defdex[:,0:5,:] = np.nan
	defdey[:,0:5,:] = np.nan
	defabs[:,0:5,:] = np.nan

	cm = f._get_defabs_cm()

	fig = plt.figure(figsize=(14.4,10.8), dpi=72)

	def update(t):
		plt.clf()

		m = Basemap(projection='npstere',boundinglat=30,lon_0=-50,resolution='l')
		x,y = m(lon,lat)
		ut,vt,xt,yt = m.transform_vector(defdex[t,::-1,:],defdey[t,::-1,:],lon[0,:],lat[::-1,0],60,60,returnxy=True)
		m.drawcoastlines()
		m.contour(x,y, oro, range(5000,80001,5000), colors='k', alpha=0.4, zorder=2)
		pabs = m.contourf(x, y, defabs[t,:,:]*86400, np.arange(4,33), zorder=1, cmap=cm)
		pang = m.quiver(xt, yt, ut, vt, zorder=3)
		panh = m.quiver(xt, yt,-ut,-vt, zorder=3)
		m.drawparallels(range(15,80,5))
		m.drawmeridians(range(0,360,30))
		plt.title(str(date_start+t*td(0.25)))
		plt.colorbar()

		return pabs.ax,

	ani = animation.FuncAnimation(fig, update, interval=500, frames=range(defabs.shape[0]))
	ani.save('test.mp4', clear_temp=False, frame_prefix=img_prefix)
	plt.show()

	return


# the end
