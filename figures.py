#!/usr/bin/python
# -*- encoding: utf-8

import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, cm as bmcm
from metopen import metopen
from utils import concat1

import static as c


# globally useful
f, oro = metopen('static', 'oro', cut=c.std_slice[1:])
oro = concat1(oro)
s   = oro.shape
lat = np.tile(f['lat'][c.std_slice[1]], (s[1],1)).T
lon = np.tile(f['lon'][c.std_slice[2]], (s[0],1))
lon = concat1(lon)
lon[:,-1] += 360
f.close()

orolevs = range(-19000,51000,2000)

# TODO: Generalisation in data fetcher <-> plotter to avoid code duplication


# #############################################################################
# 1. Plots of (multi-)yearly means for constant plev, yidx or xidx
#


# contour map of 32 year mean deformation
def map_mean_Q(q='defabs', year=None, plev=700, agg=False, quiet=False, cmap=None):
	meanZ = np.zeros(s)
	mean  = np.zeros(s)

	if year: 
		if type(year) == list:
			years = year
		elif type(year) == int:
			years = [year,]
		else:
			raise TypeError, 'year must be an integer or a list of integers'
	else:
		years = globals()['c'].years
	
	for y in years:
		if agg:	
			# TODO: Correct calculation of mean stddev!
			f, dat = metopen(c.file_stat % (y, plev, q), agg, cut=c.std_slice[1:])
		else:
			f, dat = metopen(c.file_stat % (y, plev, q), 'mean', cut=c.std_slice[1:])
		mean = concat1(dat)

		fZ, daZ = metopen(c.file_stat % (y, plev, 'Z'), 'mean', cut=c.std_slice[1:])
		meanZ += concat1(daZ)
		
		del dat, daZ
		f.close()
		fZ.close()
	
	mean  /= len(years)
	meanZ /= len(years)
	
	# Mask out surface parts below orography
	mean[meanZ < oro] = np.nan
	mean[0:5,:] = np.nan

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	#m.contourf(x,y,oro,orolevs, cmap=plt.cm.gist_earth)
	m.contourf(x, y, mean, 20, cmap=cmap)
	m.drawcoastlines()
	m.drawparallels(range(0,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# contour map of 32 year mean deformation
def wmap_mean_Q(q='defabs', year=None, plev=700, agg=False, quiet=False, cmap=None):
	meanZ = np.zeros((s[0],s[1]-1))
	mean  = np.zeros((s[0],s[1]-1))

	if year: 
		if type(year) == list:
			years = year
		elif type(year) == int:
			years = [year,]
		else:
			raise TypeError, 'year must be an integer or a list of integers'
	else:
		years = globals()['c'].years
	
	for y in years:
		if agg:	
			# TODO: Correct calculation of mean stddev!
			f, dat = metopen(c.file_stat % (y, plev, q), agg, cut=c.std_slice[1:])
		else:
			f, dat = metopen(c.file_stat % (y, plev, q), 'mean', cut=c.std_slice[1:])
		#mean = concat1(dat)
		mean += dat[::]

		fZ, daZ = metopen(c.file_stat % (y, plev, 'Z'), 'mean', cut=c.std_slice[1:])
		meanZ += daZ[::]
		
		del dat, daZ
		f.close()
		fZ.close()
	
	mean  /= len(years)
	meanZ /= len(years)
	
	# Mask out surface parts below orography
	mean[meanZ < oro[:,:-1]] = np.nan
	mean[:5,:] = np.nan
	mean[-6:,:] = np.nan

	m = Basemap(projection='robin',lon_0=0,resolution='c')
	x,y = m(lon[:,:-1],lat[:,:-1])
	#m.contourf(x,y,oro,orolevs, cmap=plt.cm.gist_earth)
	m.contourf(x, y, mean, 30, cmap=cmap)
	m.drawcoastlines()
	m.drawparallels(range(-80,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Contour map of averaged wind on top of a oro map.
def map_mean_barb(year=None, plev=700, quiver=False):
	if not year:
		fu, meanu = metopen(c.file_mstat % (plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_mstat % (plev, 'v'), 'mean', cut=c.std_slice[1:])
		fZ, meanZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
	else: 
		fu, meanu = metopen(c.file_stat % (year, plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_stat % (year, plev, 'v'), 'mean', cut=c.std_slice[1:])
		fZ, meanZ = metopen(c.file_stat % (year, plev, 'Z'), 'mean', cut=c.std_slice[1:])

	meanu[meanZ < oro[:,:-1]] = np.nan
	meanv[meanZ < oro[:,:-1]] = np.nan
	meanu = concat1(meanu[::-1,:])
	meanv = concat1(meanv[::-1,:])

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(meanu,meanv,lon[0,:],lat[::-1,0],80,80,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth)
	if not quiver:
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5)
	else:
		m.quiver(xt, yt, ut, vt)
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Contour map of averaged wind on top of a oro map.
def wmap_mean_barb(year=None, plev=700, quiver=False):
	if not year:
		fu, meanu = metopen(c.file_mstat % (plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_mstat % (plev, 'v'), 'mean', cut=c.std_slice[1:])
		fZ, meanZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
	else: 
		fu, meanu = metopen(c.file_stat % (year, plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_stat % (year, plev, 'v'), 'mean', cut=c.std_slice[1:])
		fZ, meanZ = metopen(c.file_stat % (year, plev, 'Z'), 'mean', cut=c.std_slice[1:])

	meanu[meanZ < oro[:,:-1]] = np.nan
	meanv[meanZ < oro[:,:-1]] = np.nan

	m = Basemap(projection='robin',lon_0=0,resolution='c')
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(meanu[::-1,:],meanv[::-1,:],lon[0,:-1],lat[::-1,0],80,80,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth)
	if not quiver:
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5)
	else:
		m.quiver(xt, yt, ut, vt)
	m.drawparallels(range(-80,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Contour map of averaged deformation vector on top of a oro map.
def wmap_mean_deform(year=None, plev=700):
	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mfv', cut=c.std_slice[1:])
		fZ, meanZ     = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mfv', cut=c.std_slice[1:])
		fZ, meanZ     = metopen(c.file_stat % (year, plev, 'Z'), 'mean', cut=c.std_slice[1:])
	
	meandex = np.cos(meanang[:,:]) *meanabs
	meandey = np.sin(meanang[:,:]) *meanabs
	
	meanabs[meanZ < oro[:,:-1]] = np.nan
	meanZ = meanZ[::-1,:]
	meandex[meanZ < oro[:,:-1]] = np.nan
	meandey[meanZ < oro[:,:-1]] = np.nan
	meandex[:5,:] = np.nan
	meandey[:5,:] = np.nan
	meanabs[:5,:] = np.nan
	meandex[-6:,:] = np.nan
	meandey[-6:,:] = np.nan
	meanabs[-6:,:] = np.nan

	m = Basemap(projection='robin',lon_0=0,resolution='c')
	x,y = m(lon[:,:-1],lat[:,:-1])
	ut,vt,xt,yt = m.transform_vector(meandex[::-1,:],meandey[::-1,:],lon[0,:-1],lat[::-1,0],60,60,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro[:,:-1], orolevs, cmap=plt.cm.gist_earth, zorder=1)
	m.contour(x, y, meanabs, 25, zorder=2)
	m.quiver(xt, yt, ut, vt, zorder=3)
	m.quiver(xt, yt,-ut,-vt, zorder=3)
	m.drawparallels(range(-80,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return

# Contour map of averaged deformation vector on top of a oro map.
def map_mean_deform(year=None, plev=700):
	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mfv', cut=c.std_slice[1:])
		fZ, meanZ     = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mfv', cut=c.std_slice[1:])
		fZ, meanZ     = metopen(c.file_stat % (year, plev, 'Z'), 'mean', cut=c.std_slice[1:])
	
	meandex = np.cos(meanang[:,:]) *meanabs
	meandey = np.sin(meanang[:,:]) *meanabs
	
	meanabs[meanZ < oro[:,:-1]] = np.nan
	meandex[meanZ < oro[:,:-1]] = np.nan
	meandey[meanZ < oro[:,:-1]] = np.nan
	meandex[0:5,:] = np.nan
	meandey[0:5,:] = np.nan
	meanabs[0:5,:] = np.nan

	meandex = concat1(meandex)
	meandey = concat1(meandey) 
	meanabs = concat1(meanabs)

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(meandex[::-1,:],meandey[::-1,:],lon[0,:],lat[::-1,0],60,60,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth, zorder=1)
	m.contour(x, y, meanabs, 25, zorder=2)
	m.quiver(xt, yt, ut, vt, zorder=3)
	m.quiver(xt, yt,-ut,-vt, zorder=3)
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# vertical profiles of 32years mean deformation
def ysect_mean_Q(q='defabs', year=None, yidx=57, agg=False, quiet=False, cmap=None):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	
	i = 0
	qm = np.zeros((len(c.plevs),s[1]))
	for plev in c.plevs:
		if agg:
			if year:
				f, dat = metopen(c.file_stat % (year, plev, q), agg, cut=c.std_slice[1:])
			else:
				f, dat = metopen(c.file_mstat % (plev, q), agg, cut=c.std_slice[1:])
		else:
			if year:
				f, dat = metopen(c.file_stat % (year, plev, q), 'mean', cut=c.std_slice[1:])
			else:
				f, dat = metopen(c.file_mstat % (plev, q), 'mean', cut=c.std_slice[1:])

		qm[i] = concat1(dat[yidx,:])
		i += 1
	
	orop = 1000.0*np.exp(globals()['oro'][yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	orop[orop > 1000.0] = 1000.0
	
	qm[(orop-10) < c.plevs[:,np.newaxis]] = np.nan

	plt.contour(lon[0,:], c.plevs, qm, 15, cmap=cmap)
	#plt.fill_between(lon[0,:], 1000.0-oro.values[yidx,:]/80.0, np.ones((240,))*1000.0, 'k')
	plt.fill(lon[0,:], orop, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


# 
def xsect_mean_Q(q='defabs', year=None, xidx=264, agg=False, cmap=None):
	if not quiet:
		print 'Lon %f' % lon[0,xidx]
	
	i = 0
	qm = np.zeros((len(c.plevs),s[0]))
	for plev in c.plevs:
		if agg:
			if year:
				f, dat = metopen(c.file_stat % (year, plev, q), agg)
			else:
				f, dat = metopen(c.file_mstat % (plev, q), agg)
		else:
			if year:
				f, dat = metopen(c.file_stat % (year, plev, q), 'mean')
			else:
				f, dat = metopen(c.file_mstat % (plev, q), 'mean')

		qm[i] = dat[:,xidx]
		i += 1
	
	orop = 1000.0*np.exp(globals()['oro'][:,xidx]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	orop[orop > 1000.0] = 1000.0
	
	qm[(orop-10) < c.plevs[:,np.newaxis]] = np.nan
	
	plt.contour(lat[1:-1,0], levels, qm, 15, cmap=cmap)
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return



# #############################################################################
# 2. Plots instantaneous fields or short-term averages for plev, yidx, xidx
#


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def map_date_Q(date, q='defabs', plev=700, cmap=None):
	dat = _get_instantaneous(q, date, plevs=plev)
	daZ = _get_instantaneous('Z', date, plevs=plev)

	dat[daZ < oro] = np.nan

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth)
	m.contour(x, y, dat, 25, cmap=cmap)
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Contour map of averaged deformation vector on top of a oro map.
def map_date_deform(date, plev=700):
	defabs = _get_instantaneous('defabs', date, plevs=plev)
	defang = _get_instantaneous('defang', date, plevs=plev)
	daZ    = _get_instantaneous('Z', date, plevs=plev)

	#defabs = defabs[::-1,:]
	#defang = defang[::-1,:]
	daZ = daZ[::-1,:]
	
	defdex = np.cos(defang[:,:]) *defabs
	defdey = np.sin(defang[:,:]) *defabs
	
	defabs[daZ < oro] = np.ma.masked
	defdex[daZ < oro] = np.ma.masked
	defdey[daZ < oro] = np.ma.masked
	defdex[0:5,:] = np.ma.masked
	defdey[0:5,:] = np.ma.masked
	defabs[0:5,:] = np.ma.masked

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(defdex[::-1,:],defdey[::-1,:],lon[0,:],lat[::-1,0],60,60,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth, zorder=1)
	m.contour(x, y, defabs, 25, zorder=2)
	m.quiver(xt, yt, ut, vt, zorder=3)
	m.quiver(xt, yt,-ut,-vt, zorder=3)
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Contour map of averaged wind on top of a oro map.
def map_date_barb(date, plev=700, quiver=False):
	dau = _get_instantaneous('u', date, plevs=plev)
	dav = _get_instantaneous('v', date, plevs=plev)
	daZ = _get_instantaneous('Z', date, plevs=plev)

	dau[daZ < oro] = np.nan
	dav[daZ < oro] = np.nan

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(dau[::-1,:],dav[::-1,:],lon[0,:],lat[::-1,0],80,80,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth)
	if not quiver:
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5)
	else:
		m.quiver(xt, yt, ut, vt)
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return



# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def ysect_date_Q(date, q='defabs', yidx=57, quiet=False, cmap=None):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	
	qm = _get_instantaneous(q, date, yidx=yidx, quiet=quiet)
	
	orop = 1000.0*np.exp(oro[yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	orop[orop > 1000.0] = 1000.0

	qm[orop[np.newaxis,:]-10 < c.plevs[:,np.newaxis]] = np.nan

	plt.contour(lon[0,:], c.plevs, qm, 15, cmap=cmap)
	plt.fill(lon[0,:], orop, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return



# #############################################################################
# 3. Other functions
# 

# Constant yidx and plev; plotting (multi-)year means
def ypline_mean_Q(q='defabs', yidx=57, plev=700, summarize=False, agg=False, quiet=False):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	means = {}
	for y in c.years:
		if not quiet:
			print y
		if agg:
			f, dat = metopen(c.file_stat % (y, plev, q), agg)
			means[y] = dat[yidx,:]
		else:
			f, dat = metopen(c.file_stat % (y, plev, q), 'mean')
			means[y] = dat[yidx,:]
		f.close()
	del dat
	
	if summarize:
		num = 5
		means2 = {}
		for y, dat in means.items():
			if (y-c.years[0]) % num == 0:
				key = '%d-%d' % (y,y+num-1)
				means2[key] = dat/num
			else:
				means2[key] += dat/num
		means = means2

	for y, dat in means.items():
		plt.plot(lon[0,:-1], dat, hold=True)
	
	plt.legend(means.keys())
	plt.show()

	return


# era interim orographical map 30N-90N
def map_oro():
	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	m.drawcoastlines()
	x,y = m(lon,lat)
	m.contourf(x,y, oro, orolevs, cmap=plt.cm.gist_earth)
	m.drawparallels(range(0,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return

def map_oro_dat(dat, plev=None):
	if plev:
		raise NotImplementedError
		meanabs[daZ < oro[:,:-1]] = np.nan
	
	dat = concat1(dat)
	
	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	m.drawcoastlines()
	x,y = m(lon,lat)
	m.contourf(x,y, oro, orolevs, cmap=plt.cm.gist_earth)
	m.contour(x, y, dat, 25)
	m.drawparallels(range(0,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Hovmöller diagram
def ypline_hov_Q(year, q='defabs', plev=700, yidx=57, quiet=False, cmap=None):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	f, dat = metopen(c.file_std % (year, plev, q), c.q[q])
	dat = dat[:,yidx,:]
	if f: f.close()

	dat = np.concatenate((dat[:,-60:], dat, dat[:,:60]), axis=1)
	exlon = map(lambda x: x*0.5 - 210.0, range(dat.shape[1]))
	tidxs = map(lambda x: x*0.25,range(dat.shape[0]))

	plt.contourf(exlon, tidxs, dat[:,:], 20, cmap=cmap)
	plt.plot([-180, -180], [tidxs[0], tidxs[-1]], 'k--')
	plt.plot([ 180,  180], [tidxs[0], tidxs[-1]], 'k--')
	plt.xticks(np.arange(-210,211,30))
	plt.yticks(np.array([0,31,59,90,120,151,181,212,243,273,304,334]), 
		('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
	plt.xlabel('Longitude')
	plt.ylabel(str(year))
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


def hist(year, q='defang', plev=700, yidx=57, xidx=264, quiet=False):
	if not quiet:
		print 'Lat: %f, Lon: %f' % (lat[yidx,0],lon[0,xidx])
	f, dat = metopen(c.file_stat % (year, plev, q), 'hist')
	dat = dat[:,yidx,xidx]
	if f: f.close()

	if q in c.bins:
		bins = copy.copy(c.bins[q])

		# periodic binning for defang
		if q == 'defang':
			bins[0] -= np.pi
			# TODO: Remove temporary workaround for bug fixed in dynlib.stat.binned
			dat[0], dat[-1] = dat[-1], dat[0]

		plt.bar(bins[:-1], dat[:-1], bins[1:]-bins[:-1])
		print '%d values out of binning range (%f,%f)' % (dat[-1], bins[0], bins[-1])
		plt.xlim(bins[0], bins[-1])
	else:
		raise NotImplementedError

	plt.show()

	return



# #############################################################################
# 4. Generalisations
# 


# Generalised data fetcher for instantaneous or short-term averaged fields
def _get_instantaneous(q, dates, plevs=None, yidx=None, xidx=None, tavg=True, quiet=False):
	# None means "take everything there is"
	if not plevs:
		plevs = c.plevs
	else:
		plevs = [plevs,]
	
	if yidx == None:
		yidxs = slice(None)
		ylen  = s[0]
	else:
		yidxs = yidx
		ylen  = 1
	
	if xidx == None:
		xidxs = slice(None)
		xlen  = s[1]
	else:
		xidxs = xidx
		xlen  = 1

	# Convert dates to time indexes
	if type(dates) not in ([np.ndarray, list, tuple, set]):
		dates = [dates, ]
	tidxs = map(lambda x: (x.timetuple().tm_yday-1)*4 + int(x.hour/6), dates)

	# Construct the slice
	cut = (slice(min(tidxs),max(tidxs)+1), yidxs, xidxs)
	
	# One ore more vertical levels?
	if len(plevs) > 1:
		i = 0
		dat = np.zeros((1+max(tidxs)-min(tidxs), len(c.plevs), ylen, xlen))
		#dat = dat.squeeze()
		for plev in plevs:
			if not quiet:
				print "Reading from "+c.file_std % (dates[0].year, plev, q)
			f, d = metopen(c.file_std % (dates[0].year, plev, q), c.q[q], cut=cut)
			if xidx == None:
				dat[:,i,::] = concat1(d)
			else:
				dat[:,i,::] = d
			i += 1
	else:
		if not quiet:
			print "Reading from "+c.file_std % (dates[0].year, plevs[0], q)
		f, dat = metopen(c.file_std % (dates[0].year, plevs[0], q), c.q[q], cut=cut)
		if xidx == None:
			dat = concat1(dat)
	
	# Time-averaging if specified
	if tavg and len(dates) > 1:
		dat = dat.mean(axis=0)
	
	dat = dat.squeeze()
	
	return dat


# Get aggregated (average, standard deviation, etc.) fields
def _get_aggregate(q, year=None, plev=None, yidx=None, xidx=None):


	return dat


def _get_periodic_cm():
	cdict = {'red':   ((0.0, 0.0, 0.0), (0.25, 0.8, 0.8), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0)),
		 'green': ((0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.9, 0.9), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 0.6, 0.6), (0.25, 0.0, 0.0), (0.5, 0.2, 0.2), (0.75, 0.0, 0.0), (1.0, 0.6, 0.6))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic',cdict,256)


def _get_periodic_cm2():
	cdict = {'red':   ((0.0, 0.3, 0.3), (0.25, 0.4, 0.4), (0.5, 0.6, 0.6), (0.75, 0.8, 0.8), (1.0, 0.3, 0.3)),
		 'green': ((0.0, 0.2, 0.2), (0.25, 0.8, 0.8), (0.5, 0.6, 0.6), (0.75, 0.4, 0.4), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 0.5, 0.5), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.5, 0.5))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic',cdict,256)


# that's it
