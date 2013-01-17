#!/usr/bin/env python
# -*- encoding: utf-8

import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, cm as bmcm
from datetime import datetime as dt, timedelta as td
from metopen import metopen
from utils import concat1, igauss
import windrose as wr
import streamplot as sp

import static as c
import stats

from settings import s as sts


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

scale_oro = range(10000,80001,10000)
scale_defabs = np.arange(3.0,30.1)
scale_meandefabs = np.arange(0.0,10.1,0.5)
scale_defang = (np.arange(-18,19)-0.5)*np.pi/36.0
scale_defang_coarse = np.arange(-4,5)*np.pi/8.0 - np.pi/72.0
scale_Zdiff = np.arange(-5250,5251,500)
scale_udiff = np.arange(-30,31,5)
scale_u     = np.arange(20,71,10)

ticks_defang = np.arange(-4,5)*3.1415926535/8.0 
ticklabels_defang = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']

# TODO: Generalisation in data fetcher <-> plotter to avoid code duplication


# #############################################################################
# 0. Some typically used projections
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

def NAmap():
	return Basemap(projection='lcc', lat_0=55, lat_ts=55, lon_0=-30, resolution='l', 
			width=9000000, height=6000000)
def NPmap():
	return Basemap(projection='lcc', lat_0=50, lat_ts=50, lon_0=-180, resolution='l', 
			width=9000000, height=6000000)
def Sibirmap():
	return Basemap(projection='lcc', lat_0=55, lat_ts=55, lon_0=75, resolution='l', 
			width=9000000, height=6000000)
def Ausmap():
	return Basemap(projection='lcc', lat_0=-50, lat_ts=-50, lon_0=120, resolution='l', 
			width=9000000, height=6000000)
def Bagmap():
	return Basemap(projection='lcc', lat_0=35, lat_ts=35, lon_0=55, resolution='l', 
			width=9000000, height=6000000)




# #############################################################################
# 1. Plots of (multi-)yearly means for constant plev, yidx or xidx
#


# contour map of 32 year mean deformation
def map_mean_Q(q='defabs', year=None, plev=800, agg='mean', quiet=False, cmap=None):
	if year: 
		mean  = np.zeros((s[0],s[1]-1))
		if type(year) == list:
			years = year
		elif type(year) == int:
			years = [year,]
		else:
			raise TypeError, 'year must be an integer or a list of integers'
		for y in years:
			f, dat = metopen(c.file_stat % (y, plev, q), agg, cut=c.std_slice[1:])
			mean += dat[::]
			f.close()

		del dat
		mean  /= len(years)
	else:
		f, mean = metopen(c.file_mstat % (plev, q), agg, cut=c.std_slice[1:])

	map_oro_dat(npmap(), mean, plev=plev, cmap=cmap)

	return


# contour map of trends
def map_trend_Q(q='defabs', sig=0.95, plev=800, quiet=False, cmap=None, scale=25, disable_cb=False):
	f, trend = metopen(c.file_mstat % (plev, q), 'trend', cut=c.std_slice[1:])
	ctrend = f['strend']
	overlay = map_overlay_dat(trend, scale=[0.0,99.0], colors='k', labels=False)

	upper_perc = (1.0 - sig)/2.0 + sig
	lower_perc = (1.0 - sig)/2.0
	upper_trend = trend[:,:] + ctrend[:,:]*igauss(upper_perc)
	lower_trend = trend[:,:] + ctrend[:,:]*igauss(lower_perc)

	mask = np.logical_and(upper_trend >= 0.0, lower_trend <= 0.0)
	trend[mask] = np.nan
	
	trend *= 1461.0
	if q == 'defabs':
		trend *= 1e5 #86400.0
	
	map_oro_dat(npmap(), trend, scale=scale, plev=plev, cmap=cmap, overlays=[overlay, ], disable_cb=disable_cb)

	return


# contour map of 32 year mean deformation
def wmap_mean_Q(q='defabs', year=None, plev=800, agg='mean', quiet=False, cmap=None):
	if year: 
		mean  = np.zeros((s[0],s[1]-1))
		if type(year) == list:
			years = year
		elif type(year) == int:
			years = [year,]
		else:
			raise TypeError, 'year must be an integer or a list of integers'
		for y in years:
			f, dat = metopen(c.file_stat % (y, plev, q), agg, cut=c.std_slice[1:])
			mean += dat[::]
			f.close()

		del dat
		mean /= len(years)
	else:
		f, mean = metopen(c.file_mstat % (plev, q), agg, cut=c.std_slice[1:])
	
	map_oro_dat(wmap(), mean, plev=plev, cmap=cmap)

	return


# Contour map of averaged wind on top of a oro map.
def map_mean_barb(year=None, plev=800, quiver=False):
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
def wmap_mean_barb(year=None, plev=800, quiver=False):
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
def wmap_mean_deform(year=None, plev=800):
	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	
	map_oro_deform(wmap(), meanabs, meanang, plev=plev, scale=scale_meandefabs)

	return

# Contour map of averaged deformation vector on top of a oro map.
def map_mean_deform(year=None, plev=800):
	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	
	map_oro_deform(npmap(), meanabs, meanang, plev=plev, scale=scale_meandefabs)

	return


# vertical profiles of 32years mean deformation
def ysect_mean_Q(q='defabs', year=None, yidx=51, agg=False, quiet=False, cmap=None):
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
def xsect_mean_Q(q='defabs', year=None, xidx=278, agg=False, cmap=None):
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
def map_date_Q(date, q='defabs', plev=800, cmap=None):
	dat = _get_instantaneous(q, date, plevs=plev)
	
	map_oro_dat(npmap(), dat, plev=plev, cmap=cmap)

	return


# Contour map of averaged deformation vector on top of a oro map.
def map_date_deform(date, plev=800):
	defabs = _get_instantaneous('defabs', date, plevs=plev)
	defang = _get_instantaneous('defang', date, plevs=plev)

	map_oro_deform(npmap(), defabs, defang, plev=plev)

	return


# Contour map of averaged wind on top of a oro map.
def map_date_barb(date, plev=800, quiver=False):
	dau = _get_instantaneous('u', date, plevs=plev)
	dav = _get_instantaneous('v', date, plevs=plev)
	
	map_oro_barb(npmap(), dau, dav, oro[:,:-1], 
		plev=plev, quiver=quiver, cmap=plt.cm.gist_earth, scale=orolevs)

	return


# Streamlines on a oro map.
def map_date_stream(date, plev=800, cmap=None):
	u = concat1(_get_instantaneous('u', date, plevs=plev))
	v = concat1(_get_instantaneous('v', date, plevs=plev))

	ff = np.sqrt(u*u + v*v)
	
	m = npmap()
	sp.streamplot(lon[0,:], lat[:,0], u, v, m=m, cmap=cmap)
	map_oro_dat(m, oro[:,:-1], plev=plev, cmap=plt.cm.gist_earth, scale=orolevs)

	return


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def ysect_date_Q(date, q='defabs', yidx=51, quiet=False, cmap=None):
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
def ypline_mean_Q(q='defabs', yidx=51, plev=800, summarize=False, agg=False, quiet=False):
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
def map_oro(m=npmap()):
	map_oro_dat(m, oro[:,:-1], cmap=plt.cm.gist_earth, scale=orolevs)

	return


# Hovmöller diagram
def ypline_hov_Q(year, q='defabs', plev=800, yidx=51, quiet=False, cmap=None, scale=25, slc=slice(None), 
		disable_cb=False, save=None, show=True):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	f, dat = metopen(c.file_std % (year, plev, q), c.q[q])
	fZ, daZ = metopen(c.file_std % (year, plev, 'Z'), c.q['Z'])
	dat = dat[:,yidx,:]
	daZ = daZ[:,yidx,:]
	if f: f.close()
	if fZ: fZ.close()

	oroy = oro[yidx,:-1]
	mask = oroy[np.newaxis,:] > daZ[:,:]

	if q == 'defabs':
		dat[:,:] *= 1e5 #86400.0

	dat = np.concatenate((dat[:,-60:], dat, dat[:,:60]), axis=1)
	mask = np.concatenate((mask[:,-60:], mask, mask[:,:60]), axis=1)
	dat[mask] = np.nan
	exlon = map(lambda x: x*0.5 - 210.0, range(dat.shape[1]))
	tidxs = map(lambda x: x*0.25,range(dat.shape[0]))
	tidxs = tidxs[slc]
	mask  = mask[slc,:]
	dat   = dat[slc,:]

	plt.contourf(exlon, tidxs, mask[:,:], cmap=_get_grey_cm())
	plt.contourf(exlon, tidxs, dat[:,:], scale, colors=_sample_cm(cmap, len(scale)-1), extend='both')
	plt.plot([-180, -180], [tidxs[0], tidxs[-1]], 'k--')
	plt.plot([ 180,  180], [tidxs[0], tidxs[-1]], 'k--')
	plt.xticks(np.arange(-210,211,30))
	plt.yticks(np.array([0,31,59,90,120,151,181,212,243,273,304,334]), 
		('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
	plt.xlabel('Longitude')
	plt.ylabel(str(year))
	if slc.start:
		plt.ylim((slc.start/4.0, plt.ylim()[1]))
	if slc.stop:
		plt.ylim((plt.ylim()[0], slc.stop/4.0))
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	if not disable_cb:
		plt.colorbar(pad=0.04, fraction=0.04)
	
	
	if save:
		plt.savefig(save, format='png')
	if show:
		plt.show()

	return


def hist(year, q='defang', plev=800, yidx=51, xidx=278, quiet=False):
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
			#dat[0], dat[-1] = dat[-1], dat[0]

		plt.bar(bins[:-1], dat[:-1], bins[1:]-bins[:-1], zorder=1)
		plt.plot((bins[:-1]+bins[1:])/2.0, stats.running_mean_periodic(dat[:-1], rlen=7), 'r', zorder=2)
		print '%d values out of binning range (%f,%f)' % (dat[-1], bins[0], bins[-1])
		plt.xlim(bins[0], bins[-1])
	else:
		raise NotImplementedError

	plt.show()

	return


def phist(year, q='defang', plev=800, yidx=51, xidx=278, quiet=False):
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
			#dat[0], dat[-1] = dat[-1], dat[0]

		print '%d values out of binning range (%f,%f)' % (dat[-1], bins[0], bins[-1])
		wr.bar(dat[:-1], np.array([0,]), bins=1, prebinned=bins, halfplot=True, normed=True)
		#wr.contour(stats.running_mean_periodic(dat[:-1], rlen=7), 
		#		np.array([0,]), bins=1, prebinned=bins, halfplot=True, normed=True, colors='r')

	else:
		raise NotImplementedError

	plt.show()

	return


def ts_phist(qd='defang', qv='defabs', plev='800', pos='Greenland_TB', mons=[], years=[], mask=None,
		quiet=False, show=True, save='', title='', cmap=None, disable_legend=False):
	fd, dd = metopen('../timeseries/%s.%s.%s_ts' % (pos, plev, qd), 'ts', cut=(slice(None),) )
	fv, dv = metopen('../timeseries/%s.%s.%s_ts' % (pos, plev, qv), 'ts', cut=(slice(None),) )
	if not quiet:
		yidx, xidx = fd['pos']
		print '%s (Lat: %f, Lon: %f)' % (pos, lat[yidx,0],lon[0,xidx])
	fd.close()
	fv.close()
	
	if (len(mons) > 0 or len(years) > 0) and not mask:
		dates = [dt(1979,1,1,0)+td(0.25)*i for i in range(len(dd))]

	if type(mask) == np.ndarray:
		pass
	elif len(mons) > 0 and len(years) > 0:
		mask = np.array(map(lambda x: x.month in mons and x.year in years, dates))
	elif len(mons) > 0:
		mask =  np.array(map(lambda x: x.month in mons, dates))
	elif len(years) > 0:
		mask =  np.array(map(lambda x: x.year in years, dates))
	
	if type(mask) == np.ndarray:
		dd = dd[mask]
		dv = dv[mask]
	
	print len(dd)
	
	# wr is made for wind direction: N -> 0, clockwise -> pos; 
	# however for deformation angle: E -> 0, anticlockwise -> pos;
	dd = (90.0 - dd*(180/np.pi)) % 360.0
	dv *= 1e5 #86400.0
	
	# doubling
	dd = np.append(dd, (dd + 180.0) % 360.0)
	dv = np.append(dv, dv)
	
	wr.bar(dd, dv, bins=[0, 6, 12, 18], halfplot=True, normed=True, nsector=72, 
			disable_legend=disable_legend, cmap=cmap)
	plt.text(1.65*np.pi, 1.2*plt.gca().get_rmax(), 'Pos: %3.1f%s, %3.1f%s' % (
			abs(lat[yidx,0]), 'N' if lat[yidx,0] > 0 else 'S', 
			abs(lon[0,xidx]), 'E' if lon[0,xidx] > 0 else 'W'),
		fontsize=14)
	
	if title:
		plt.title(title)
	if save:
		plt.savefig(save, format='png')
	if show:
		plt.show()

	return


def ts_hist(q='defabs', plev='800', pos='Greenland_TB', bins=20, mons=[], years=[], mask=None, ylim=None,
		quiet=False, show=True, save='', title='', color=None):
	f, dat = metopen('../timeseries/%s.%s.%s_ts' % (pos, plev, q), 'ts', cut=(slice(None),) )
	if not quiet:
		yidx, xidx = f['pos']
		print '%s (Lat: %f, Lon: %f)' % (pos, lat[yidx,0],lon[0,xidx])
	f.close()
	
	if (len(mons) > 0 or len(years) > 0) and not mask:
		dates = [dt(1979,1,1,0)+td(0.25)*i for i in range(len(dat))]

	if type(mask) == np.ndarray:
		pass
	elif len(mons) > 0 and len(years) > 0:
		mask = np.array(map(lambda x: x.month in mons and x.year in years, dates))
	elif len(mons) > 0:
		mask =  np.array(map(lambda x: x.month in mons, dates))
	elif len(years) > 0:
		mask =  np.array(map(lambda x: x.year in years, dates))
	
	if type(mask) == np.ndarray:
		dat = dat[mask]
	
	print len(dat)
	
	if q == 'defabs':
		dat *= 1e5 #86400.0
		if not ylim: ylim = 0.3
	
	plt.hist(dat, bins=bins, normed=True, color=color)
	plt.annotate('Pos: %3.1f%s, %3.1f%s' % (
			abs(lat[yidx,0]), 'N' if lat[yidx,0] > 0 else 'S', 
			abs(lon[0,xidx]), 'E' if lon[0,xidx] > 0 else 'W'),
		xy=(0.7,0.92), fontsize=14, xycoords='axes fraction')

	if ylim:
		plt.ylim(0,ylim)
	
	if title:
		plt.title(title)
	if save:
		plt.savefig(save, format='png')
	if show:
		plt.show()

	return


def ts_wavelet(q='defabs', plev='800', pos='Greenland_TB', scale=25, cmap=None, p=np.pi,
		quiet=False, show=True, save='', title=''):
	import mlpy

	f, ts = metopen('../timeseries/%s.%s.%s_ts' % (pos, plev, q), 'ts', cut=(slice(None),) )
	if not quiet: 
		yidx, xidx = f['pos']
		print '%s (Lat: %f, Lon: %f)' % (pos, lat[yidx,0],lon[0,xidx])
	f.close()

	if q == 'defabs':
		ts *= 1e5 #86400.0

	ts_2011 = ts[-1460:]
	tswl_2011, scales_2011 = mlpy.cwt(x=ts_2011, dt=1, dj=0.125, wf='morlet', p=p)
	dates_2011 = [dt(2011,1,1,0)+td(0.25)*i for i in range(1460)]
	
	ts_weekly = np.array(map(lambda i: ts[i*28:(i+1)*28].mean(), range(len(ts)/28)))
	tswl_weekly, scales_weekly = mlpy.cwt(x=ts_weekly, dt=1, dj=0.125, wf='morlet', p=p)
	dates_weekly = [dt(1979,1,1,0)+td(7)*i for i in range(len(ts_weekly))]
	
	tswl_mean = np.zeros(tswl_2011.shape)
	dates = [dt(1979,1,1,0)+td(0.25)*i for i in range(len(ts))]
	for yr in range(1979,2012):
		tidx = dates.index(dt(yr,1,1,0))
		tswl_out, scales_mean = mlpy.cwt(x=ts[tidx:tidx+1460], dt=1, dj=0.125, wf='morlet', p=p)
		tswl_mean += np.abs(tswl_out)
	tswl_mean /= 33.0
	
	# normalisation of the scale for different shape parameters p
	norm = 2.0*np.pi/p

	fig = plt.gcf()
	dfmt = mdates.DateFormatter('%b')
	ax1 = fig.add_subplot(211)
	plt.contourf(dates_weekly, scales_weekly*7.0/365.25*norm, np.abs(tswl_weekly), scale, cmap=cmap, extend='max' )
	ax1.set_yscale('log')
	ax1.set_ylabel('Period [years]')
	plt.colorbar()
	#ax2 = fig.add_subplot(312)
	#plt.contourf(dates_2011, scales_2011*0.25*norm, np.abs(tswl_2011), scale, cmap=cmap, extend='max' )
	#ax2.set_yscale('log')
	#ax2.set_ylabel('Period [days]')
	#ax2.xaxis.set_major_formatter(dfmt)
	#plt.colorbar()
	ax3 = fig.add_subplot(212)
	plt.contourf(dates_2011, scales_mean*0.25*norm, tswl_mean, scale, cmap=cmap, extend='max' )
	ax3.set_yscale('log')
	ax3.set_ylabel('Period [days]')
	ax3.xaxis.set_major_formatter(dfmt)
	plt.colorbar()
	
	if title:
		plt.title(title)
	if save:
		plt.savefig(save, format='png')
	if show:
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
		xlen  = s[1]-1
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

def _get_greys_cm():
	cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'green': ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'blue':  ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1))  }

	return mpl.colors.LinearSegmentedColormap('my_grey',cdict,256)

def _get_defabs_cm():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.4, 0.4), (0.867, 1.0, 1.0), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.5, 0.5), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 1.0, 1.0), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def _get_defabs_cm2():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def _get_q_cm():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.30, 0.30),  (0.867, 0.1, 0.1), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.65, 0.65),  (0.867, 0.2, 0.2), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 0.80, 0.80),  (0.867, 0.6, 0.6), (1.0, 0.8, 0.8))  }

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


def _get_periodic_cm3():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.03, 0.9, 0.9), (0.27, 0.3, 0.2), (0.515, 0.5, 0.5), (0.76, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.03, 0.9, 0.9), (0.27, 0.3, 0.2), (0.515, 0.5, 0.5), (0.76, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.03, 1.0, 1.0), (0.27, 0.7, 0.7), (0.515, 0.5, 0.5), (0.76, 0.20, 0.3), (1.0, 0.9, 1.0)) }

	return mpl.colors.LinearSegmentedColormap('my_periodic3',cdict,256)



# #############################################################################
# 6. Generalised data plotters
# 


def map_oro_dat(m, dat, **kwargs):
	# TODO: q as argument or merging kwargs earlier?
	kwargs = sts.contourf.merge('defabs', **kwargs)

	dat = concat1(dat)
	plev = kwargs.pop('plev')
	if plev:
		f,daZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
		if f: f.close()
		daZ = concat1(daZ)
		mask = daZ[:,:] < oro[:,:]
		dat[mask] = np.nan
	else:
		mask = slice(None)
	
	m.drawcoastlines(color=kwargs.pop('coastcolor'))
	x,y = m(lon,lat)
	m.contour(x,y, oro, scale_oro, colors=kwargs.pop('orocolor'), alpha=kwargs.pop('oroalpha'), zorder=2)
	if plev:
		m.contourf(x, y, mask, colors=kwargs.pop('maskcolor'))
	
	scale = kwargs.pop('scale')
	cs = m.contourf(x, y, dat, scale, zorder=1, **kwargs)
	if not type(scale) == int:
		cs.set_clim(scale[0], scale[-1])

	gridcolor = kwargs.pop('gridcolor')
	m.drawparallels(range(-80,81,5), color=gridcolor)
	m.drawmeridians(range(0,360,30), color=gridcolor)

	if not kwargs.pop('disable_cb'):
		cb = plt.colorbar(ticks=kwargs.pop('ticks'), shrink=0.85, pad=0.015, fraction=0.10)
		if kwargs.get('ticklabels'): 
			cb.ax.set_yticklabels(kwargs.pop('ticklabels'))
	
	for overlay in kwargs.pop('overlays'):
		overlay(m,x,y, zorder=2, mask=mask)
	
	if kwargs.get('mark'):
		yidx, xidx = kwargs.pop('mark')
		m.scatter(x[yidx,xidx], y[yidx,xidx], 484, marker='o', facecolors=(0,0,0,0), 
				edgecolors='k', linewidths=3, zorder=3)
	
	if kwargs.get('title'):
		plt.title(kwargs.pop('title'))
	if kwargs.get('save'):
		plt.savefig(kwargs.pop('save'), format='png')
	if kwargs.pop('show'):
		plt.show()

	return


def map_oro_deform(m, defabs, defang, **kwargs):
	kwargs = sts.contourf.merge('defabs', **kwargs)

	defabs = concat1(defabs*1e5)
	defang = concat1(defang)
	defdex = np.cos(defang[:,:]) *defabs
	defdey = np.sin(defang[:,:]) *defabs
	plev = kwargs.pop('plev')
	if plev and not type(daZ) == np.ndarray:
		f,daZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
		if f: f.close()
	if type(daZ) == np.ndarray:
		daZ = concat1(daZ)
		mask = daZ[:,:] < oro[:,:]
		defabs[mask] = np.nan
		defang[mask] = np.nan
		defdex[mask] = np.nan
		defdey[mask] = np.nan
	else:
		mask = slice(None)
	
	m.drawcoastlines(color=kwargs.pop('coastcolor'))
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(defdex[::-1,:],defdey[::-1,:],lon[0,:],lat[::-1,0], 24, 16, returnxy=True)
	m.contour(x,y, oro, scale_oro, colors=kwargs.pop('orocolor'), alpha=kwargs.pop('oroalpha'), zorder=2)
	if type(daZ) == np.ndarray:
		m.contourf(x, y, mask, colors=kwargs.pop('maskcolor'))
	
	scale = kwargs.pop('scale')
	cs = m.contourf(x, y, defabs, scale, zorder=1, **kwargs)
	if not type(scale) == int:
		cs.set_clim(scale[0], scale[-1])
	m.quiver(xt, yt, ut, vt, zorder=4, scale=360, alpha=0.7)
	m.quiver(xt, yt, -ut, -vt, zorder=4, scale=360, alpha=0.7)
	
	gridcolor = kwargs.pop('gridcolor')
	m.drawparallels(range(-80,81,10), color=gridcolor)
	m.drawmeridians(range(0,360,30), color=gridcolor)

	if not kwargs.pop('disable_cb'):
		plt.colorbar(ticks=kwargs.pop('ticks'), orientation='horizontal', shrink=0.8, fraction=0.08, pad=0.02)
		if kwargs.get('ticklabels'): 
			cb.ax.set_yticklabels(kwargs.pop('ticklabels'))
	
	for overlay in kwargs.pop('overlays'):
		overlay(m,x,y, zorder=3, mask=mask)
	
	if kwargs.get('mark'):
		yidx, xidx = kwargs.pop('mark')
		m.scatter(x[yidx,xidx], y[yidx,xidx], 484, marker='o', facecolors=(0,0,0,0), 
				edgecolors='k', linewidths=3, zorder=3)
	
	if kwargs.pop('title'):
		plt.title(title)
	if kwargs.pop('save'):
		plt.savefig(save, format='png')
	if kwargs.pop('show'):
		plt.show()

	return


def map_oro_barb(m, u, v, dat=None, **kwargs):
	kwargs = sts.contourf.merge('defabs', **kwargs)

	u   = concat1(u)
	v   = concat1(v)
	if not dat == None: dat = concat1(dat)
	
	plev = kwargs.pop('plev')
	if plev:
		f,daZ = metopen(c.file_mstat % (plev, 'Z'), 'mean', cut=c.std_slice[1:])
		if f: f.close()
		daZ = concat1(daZ)
		u[daZ < oro[:,:]] = np.nan
		v[daZ < oro[:,:]] = np.nan
		if not dat == None: dat[daZ < oro[:,:]] = np.nan
	
	m.drawcoastlines(color=kwargs.pop('coastcolor'))
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(u[::-1,:],v[::-1,:],lon[0,:],lat[::-1,0], 60, 60, returnxy=True)
	m.contour(x,y, oro, scale_oro, colors=kwargs.pop('orocolor'), alpha=kwargs.pop('oroalpha'), zorder=2)

	if plev:
		m.contourf(x, y, daZ < oro[:,:], colors=kwargs.pop('maskcolor'))
	
	scale = kwargs.pop('scale')
	if not dat == None: m.contourf(x, y, dat, scale, zorder=1, **kwargs)
	if not quiver:
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5, zorder=3)
	else:
		m.quiver(xt, yt, ut, vt, zorder=3)
	
	gridcolor = kwargs.pop('gridcolor')
	m.drawparallels(range(-80,81,5), color=gridcolor)
	m.drawmeridians(range(0,360,30), color=gridcolor)

	if not kwargs.pop('disable_cb'):
		cb = plt.colorbar(ticks=kwargs.pop('ticks'), shrink=0.85, pad=0.015, fraction=0.10)
		if kwargs.get('ticklabels'): 
			cb.ax.set_yticklabels(kwargs.pop('ticklabels'))
	
	if kwargs.get('mark'):
		yidx, xidx = kwargs.pop('mark')
		m.scatter(x[yidx,xidx], y[yidx,xidx], 484, marker='o', facecolors=(0,0,0,0), 
				edgecolors='k', linewidths=3, zorder=3)
	
	if kwargs.pop('title'):
		plt.title(title)
	if kwargs.pop('save'):
		plt.savefig(save, format='png')
	if kwargs.pop('show'):
		plt.show()

	return


def map_overlay_dat(dat, **kwargs):  
	dat = concat1(dat)

	def overlay(m, x, y, zorder, mask=None):
		if type(mask) == np.ndarray:
			dat[mask] = np.nan
		cs =  m.contour(x, y, dat, **kwargs)
		#if labels:
		#	plt.clabel(cs, fontsize=12, inline=True, inline_spacing=2)

		return

	return overlay


def map_overlay_mask(dat, **kwargs):  
	dat = concat1(dat)

	def overlay(m, x, y, zorder, mask=None):
		if type(mask) == np.ndarray:
			raise NotImplementedError, 'here be dragons'
		cs =  m.contourf(x, y, dat, **kwargs)

		return

	return overlay


def phist_defang(defang, defabs=None):
	pass


# that's it
