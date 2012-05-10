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
import stats


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

wmap_scale_oro = range(10000,80001,10000)
wmap_scale_defabs = np.arange(3.0,30.1)
wmap_scale_meandefabs = np.arange(0.0,10.1,0.5)

# TODO: Generalisation in data fetcher <-> plotter to avoid code duplication


# #############################################################################
# 1. Plots of (multi-)yearly means for constant plev, yidx or xidx
#


# contour map of 32 year mean deformation
def map_mean_Q(q='defabs', year=None, plev=700, agg=False, quiet=False, cmap=None):
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
		mean = dat

		del dat
		f.close()
	
	mean  /= len(years)

	map_oro_dat(wmap(), mean, plev=plev, cmap=cmap)

	return


# contour map of 32 year mean deformation
def wmap_mean_Q(q='defabs', year=None, plev=700, agg=False, quiet=False, cmap=None):
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
		
		del dat
		f.close()
	
	mean  /= len(years)
	
	map_oro_dat(wmap(), mean, plev=plev, cmap=cmap)

	return


# Contour map of averaged wind on top of a oro map.
def map_mean_barb(year=None, plev=700, quiver=False):
	if not year:
		fu, meanu = metopen(c.file_mstat % (plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_mstat % (plev, 'v'), 'mean', cut=c.std_slice[1:])
	else: 
		fu, meanu = metopen(c.file_stat % (year, plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_stat % (year, plev, 'v'), 'mean', cut=c.std_slice[1:])

	map_oro_barb(npmap(), meanu, meanv, oro[:,:-1], 
		plev=plev, quiver=quiver, cmap=plt.cm.gist_earth, scale=orolevs)

	return


# Contour map of averaged wind on top of a oro map.
def wmap_mean_barb(year=None, plev=700, quiver=False):
	if not year:
		fu, meanu = metopen(c.file_mstat % (plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_mstat % (plev, 'v'), 'mean', cut=c.std_slice[1:])
	else: 
		fu, meanu = metopen(c.file_stat % (year, plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_stat % (year, plev, 'v'), 'mean', cut=c.std_slice[1:])

	map_oro_barb(wmap(), meanu, meanv, oro[:,:-1], 
		plev=plev, quiver=quiver, cmap=plt.cm.gist_earth, scale=orolevs)

	return


# Contour map of averaged deformation vector on top of a oro map.
def wmap_mean_deform(year=None, plev=700):
	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	
	map_oro_deform(wmap(), meanabs, meanang, plev=plev, scale=wmap_scale_meandefabs)

	return

# Contour map of averaged deformation vector on top of a oro map.
def map_mean_deform(year=None, plev=700):
	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	
	map_oro_deform(npmap(), meanabs, meanang, plev=plev, scale=wmap_scale_meandefabs)

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
	
	map_oro_dat(npmap(), dat, plev=plev, cmap=cmap)

	return


# Contour map of averaged deformation vector on top of a oro map.
def map_date_deform(date, plev=700):
	defabs = _get_instantaneous('defabs', date, plevs=plev)
	defang = _get_instantaneous('defang', date, plevs=plev)

	map_oro_deform(npmap(), defabs, defang, plev=plev)

	return


# Contour map of averaged wind on top of a oro map.
def map_date_barb(date, plev=700, quiver=False):
	dau = _get_instantaneous('u', date, plevs=plev)
	dav = _get_instantaneous('v', date, plevs=plev)
	
	map_oro_barb(npmap(), dau, dav, oro[:,:-1], 
		plev=plev, quiver=quiver, cmap=plt.cm.gist_earth, scale=orolevs)

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


# era interim orographical map 15N-90N
def map_oro():
	map_oro_dat(npmap(), oro[:,:-1], cmap=plt.cm.gist_earth, scale=orolevs)

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

		plt.bar(bins[:-1], dat[:-1], bins[1:]-bins[:-1], zorder=1)
		plt.plot((bins[:-1]+bins[1:])/2.0, stats.running_mean_periodic(dat[:-1], rlen=7), 'r', zorder=2)
		print '%d values out of binning range (%f,%f)' % (dat[-1], bins[0], bins[-1])
		plt.xlim(bins[0], bins[-1])
	else:
		raise NotImplementedError

	plt.show()

	return



# #############################################################################
# 4. Generalised data fetchers
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
			dat[:,i,::] = d
			i += 1
	else:
		if not quiet:
			print "Reading from "+c.file_std % (dates[0].year, plevs[0], q)
		f, dat = metopen(c.file_std % (dates[0].year, plevs[0], q), c.q[q], cut=cut)
	
	# Time-averaging if specified
	if tavg and len(dates) > 1:
		dat = dat.mean(axis=0)
	
	dat = dat.squeeze()
	
	return dat


# Get aggregated (average, standard deviation, etc.) fields
def _get_aggregate(q, year=None, plev=None, yidx=None, xidx=None):


	return dat


# #############################################################################
# 5. Colour maps
# 

def _get_grey_cm():
	cdict = {'red':   ((0.0, 0.4, 0.4), (1.0, 0.4, 0.4)),
		 'green': ((0.0, 0.4, 0.4), (1.0, 0.4, 0.4)),
		 'blue':  ((0.0, 0.4, 0.4), (1.0, 0.4, 0.4))  }

	return mpl.colors.LinearSegmentedColormap('my_grey',cdict,256)


def _get_defabs_cm():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.4, 0.4), (0.867, 1.0, 1.0), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.5, 0.5), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 1.0, 1.0), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)


def _get_periodic_cm():
	cdict = {'red':   ((0.0, 0.0, 0.0), (0.25, 0.8, 0.8), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0)),
		 'green': ((0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.9, 0.9), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 0.6, 0.6), (0.25, 0.0, 0.0), (0.5, 0.2, 0.2), (0.75, 0.0, 0.0), (1.0, 0.6, 0.6))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic',cdict,256)


def _get_periodic_cm2():
	cdict = {'red':   ((0.0, 0.3, 0.3), (0.25, 0.4, 0.4), (0.5, 0.6, 0.6), (0.75, 0.8, 0.8), (1.0, 0.3, 0.3)),
		 'green': ((0.0, 0.2, 0.2), (0.25, 0.8, 0.8), (0.5, 0.6, 0.6), (0.75, 0.4, 0.4), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 0.5, 0.5), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.5, 0.5))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic2',cdict,256)


# #############################################################################
# 6. Some typically used projections
# 


def wmap():
	return Basemap(projection='robin',lon_0=0,resolution='c')

def npmap():
	return Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')

def spmap():
	return Basemap(projection='spstere',boundinglat=-15,lon_0=0,resolution='l')

def Greenmap():
	return Basemap(projection='stere', lat_0=65, lat_ts=65, lon_0=-50, resolution='l', 
			width=8100000, height=5400000)

# #############################################################################
# 7. Generalised data plotters
# 


def map_oro_dat(m, dat, plev=None, mark=None, cmap=None, scale=25):
	dat = concat1(dat)
	if plev:
		f,daZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
		if f: f.close()
		daZ = concat1(daZ)
		dat[daZ < oro[:,:]] = np.nan
	
	#dat[:5,:] = np.nan
	#dat[-5:,:] = np.nan

	m.drawcoastlines()
	x,y = m(lon,lat)
	m.contour(x,y, oro, wmap_scale_oro, colors='k', alpha=0.4, zorder=2)
	if plev:
		m.contourf(x, y, daZ < oro[:,:], cmap=_get_grey_cm())
	m.contourf(x, y, dat, scale, zorder=1, cmap=cmap, extend='both')
	m.drawparallels(range(-80,81,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	if mark:
		yidx, xidx = mark
		m.scatter(x[yidx,xidx], y[yidx,xidx], 49, marker='+', color='r', zorder=3)
	plt.show()

	return


def map_oro_deform(m, defabs, defang, plev=None, mark=None, cmap=_get_defabs_cm(), scale=wmap_scale_defabs):
	defabs = concat1(defabs)
	defang = concat1(defang)
	defdex = np.cos(defang[:,:]) *defabs
	defdey = np.sin(defang[:,:]) *defabs
	if plev:
		f,daZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
		if f: f.close()
		daZ = concat1(daZ)
		defabs[daZ < oro[:,:]] = np.nan
		defang[daZ < oro[:,:]] = np.nan
		defdex[daZ < oro[:,:]] = np.nan
		defdey[daZ < oro[:,:]] = np.nan
	
	m.drawcoastlines()
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(defdex[::-1,:],defdey[::-1,:],lon[0,:],lat[::-1,0],50,50,returnxy=True)
	m.contour(x,y, oro, wmap_scale_oro, colors='k', alpha=0.4, zorder=2)
	if plev:
		m.contourf(x, y, daZ < oro[:,:], cmap=_get_grey_cm())
	m.contourf(x, y, defabs*86400, scale, zorder=1, cmap=cmap, extend='both')
	m.quiver(xt, yt, ut, vt, zorder=3)
	m.quiver(xt, yt, -ut, -vt, zorder=3)
	m.drawparallels(range(-80,81,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	if mark:
		yidx, xidx = mark
		m.scatter(x[yidx,xidx], y[yidx,xidx], 25, marker='+', color='r')
	plt.show()

	return


def map_oro_barb(m, u, v, dat=None, plev=None, mark=None, quiver=False, cmap=None, scale=25):
	u   = concat1(u)
	v   = concat1(v)
	if not dat == None: dat = concat1(dat)
	if plev:
		f,daZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
		if f: f.close()
		daZ = concat1(daZ)
		u[daZ < oro[:,:]] = np.nan
		v[daZ < oro[:,:]] = np.nan
		if not dat == None: dat[daZ < oro[:,:]] = np.nan
	
	m.drawcoastlines()
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(u[::-1,:],v[::-1,:],lon[0,:],lat[::-1,0],60,60,returnxy=True)
	m.contour(x,y, oro, wmap_scale_oro, colors='k', alpha=0.4, zorder=2)
	if plev:
		m.contourf(x, y, daZ < oro[:,:], cmap=_get_grey_cm())
	if not dat == None: m.contourf(x, y, dat, scale, cmap=cmap, zorder=1, extend='both')
	if not quiver:
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5, zorder=3)
	else:
		m.quiver(xt, yt, ut, vt, zorder=3)
	m.drawparallels(range(-80,81,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	if mark:
		yidx, xidx = mark
		m.scatter(x[yidx,xidx], y[yidx,xidx], 25, marker='+', color='r')
	plt.show()

	return


# that's it
