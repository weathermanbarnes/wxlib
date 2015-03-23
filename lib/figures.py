#!/usr/bin/env python
# -*- encoding: utf-8

import copy
import numpy as np
import scipy.interpolate as intp
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import mpl_toolkits.basemap
import mpl_toolkits.basemap.cm as bmcm
Basemap = mpl_toolkits.basemap.Basemap
basemap_version = mpl_toolkits.basemap.__version__

from datetime import datetime as dt, timedelta as td
from metopen import metopen
from utils import concat1, concat1lonlat, igauss, __unflatten_fronts_t, sect_gen_points
import windrose as wr
import streamplot as sp

import stats

from settings import conf as c
from autoscale import autoscale


# TODO: Generalisation i0n data fetcher <-> plotter to avoid code duplication
# TODO: Cleanup import!
# TODO: avoid that line!
#s = (361,721)


# #############################################################################
# 1. Plots of (multi-)yearly means for constant plev, yidx or xidx
#

# contour map of 32 year mean deformation
def map_mean_Q(q, agg='mean', year=None, **kwargs):
	kwargs = c.contourf.merge(q, **kwargs)
	plev = kwargs.get('plev')

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
	
	map_oro_dat(kwargs.pop('m'), mean, **kwargs)

	return


# Contour map of averaged wind on top of a oro map.
def map_mean_barb(q='oro', year=None, quiver=False, **kwargs):
	kwargs = c.contourf.merge(q, **kwargs)
	plev = kwargs.get('plev')

	if not year:
		fu, meanu = metopen(c.file_mstat % (plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_mstat % (plev, 'v'), 'mean', cut=c.std_slice[1:])
	else: 
		fu, meanu = metopen(c.file_stat % (year, plev, 'u'), 'mean', cut=c.std_slice[1:])
		fv, meanv = metopen(c.file_stat % (year, plev, 'v'), 'mean', cut=c.std_slice[1:])
	if q == 'oro':
		dat = oro[:,:-1]
	else:
		if not year:
			f, dat = metopen(c.file_mstat % (plev, q), 'mean', cut=c.std_slice[1:])
		else:
			f, dat = metopen(c.file_stat % (year, plev, q), 'mean', cut=c.std_slice[1:])

	kwargs['quiver'] = quiver
	map_oro_barb(kwargs.pop('m'), meanu, meanv, dat, **kwargs)

	return


# Contour map of averaged deformation vector on top of a oro map.
def map_mean_deform(year=None, **kwargs):
	kwargs = c.contourf.merge('defabs', **kwargs)
	kwargs['scale'] = settings.scale_defabs_mean
	kwargs['extend'] = 'both'
	plev = kwargs.get('plev')

	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean', cut=c.std_slice[1:])
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mfv', cut=c.std_slice[1:])
	
	map_oro_deform(kwargs.pop('m'), meanabs, meanang, **kwargs)

	return


# contour map of trends
def map_trend_Q(q, sig=0.95, **kwargs):
	kwargs = c.contourf.merge(q, **kwargs)
	plev = kwargs.get('plev')

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
	kwargs = c.contourf.merge(q, **kwargs)
	kwargs['overlay'] = [overlay, ]
	map_oro_dat(kwargs.pop('m'), trend, **kwargs)

	return


# vertical profiles of 32years mean deformation
# TODO: Generalize to accept kwargs!
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
	
	if q in settings.hooks:
		qm = settings.hooks[q](qm)

	qm[(orop-10) < c.plevs[:,np.newaxis]] = np.nan

	plt.contour(lon[0,:], c.plevs, qm, 15, cmap=cmap)
	#plt.fill_between(lon[0,:], 1000.0-oro.values[yidx,:]/80.0, np.ones((240,))*1000.0, 'k')
	plt.fill(lon[0,:], orop, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


# vertical profiles of 32years mean deformation
# TODO: Generalize to accept kwargs!
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
	
	if q in settings.hooks:
		qm = settings.hooks[q](qm)

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
def map_date_Q(q, date, **kwargs):
	dat = _get_instantaneous(q, date, kwargs.get('plev'))
	
	kwargs = c.contourf.merge(q, **kwargs)
	map_oro_dat(kwargs.pop('m'), dat, **kwargs)

	return


# Contour map of averaged deformation vector on top of a oro map.
def map_date_deform(date, **kwargs):
	defabs = _get_instantaneous('defabs', date, kwargs.get('plev'))
	defang = _get_instantaneous('defang', date, kwargs.get('plev'))

	kwargs = c.contourf.merge(q, **kwargs)
	map_oro_deform(kwargs.pop('m'), defabs, defang, **kwargs)

	return


# Contour map of averaged wind on top of a oro map.
def map_date_barb(date, q='oro', quiver=False, **kwargs):
	dau = _get_instantaneous('u', date, kwargs.get('plev'))
	dav = _get_instantaneous('v', date, kwargs.get('plev'))

	if q == 'oro':
		dat = oro[:,:-1]
	else:
		dat = _get_instantaneous(q, date, kwargs.get('plev'))
	
	kwargs = c.contourf.merge(q, **kwargs)
	kwargs['quiver'] = quiver
	map_oro_barb(kwargs.pop('m'), dau, dav, dat, **kwargs)

	return


# Streamlines on a contourf map.
def map_date_stream(date, q='oro', **kwargs):
	u = concat1(_get_instantaneous('u', date, plevs=plev))
	v = concat1(_get_instantaneous('v', date, plevs=plev))

	if q == 'oro':
		dat = oro[:,:-1]
	elif q == 'ff':
		dat = ff = np.sqrt(u*u + v*v)
	else:
		dat = _get_instantaneous(q, date, kwargs.get('plev'))
	
	if q in settings.hooks:
		dat = settings.hooks[q](dat)

	m = kwargs.pop('m')
	kwargs = c.contourf.merge(q, **kwargs)
	sp.streamplot(lon[0,:], lat[:,0], u, v, m=m, **s.contour.u)
	map_oro_dat(m, dat, **kwargs)

	return


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
# TODO: Generalize to accept kwargs!
def ysect_date_Q(date, q='defabs', yidx=51, quiet=False, cmap=None):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	
	qm = _get_instantaneous(q, date, yidx=yidx, quiet=quiet)
	
	orop = 1000.0*np.exp(oro[yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	orop[orop > 1000.0] = 1000.0

	if q in settings.hooks:
		qm = settings.hooks[q](qm)

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
# TODO: Generalize to accept kwargs!
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
def map_oro(**kwargs):
	kwargs = c.contourf.merge('oro', **kwargs)
	map_oro_dat(kwargs.pop('m'), oro[:,:-1], **kwargs)

	return


# Hovmöller diagram
# TODO: Generalize to accept kwargs!
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


# TODO: Generalize to accept kwargs!
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


# TODO: Generalize to accept kwargs!
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


# TODO: Generalize to accept kwargs!
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


# TODO: Generalize to accept kwargs!
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


# TODO: Generalize to accept kwargs!
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
# 4. Generalised data plotters
# 
def sect_oro_dat(dat, ps, sect, static, datmap=None, p=None, **kwargs):
	# 1a. Prepare
	kwargs = __prepare_config(kwargs)
	kwargs_map = copy.copy(kwargs)
	mask = __map_create_mask(static, kwargs)
	
	if 'map_axes' in kwargs:
		plt.sca(kwargs['map_axes'])
	else:
		plt.subplot(212)
	m, x, y, lon, lat = __map_setup(mask, static, kwargs)
	
	# 1b. Create and interpolate points of the cross section
	xlon, xlat, xxy = sect_gen_points(sect, m, 25000.0)

	dati = np.empty((dat.shape[0], len(xlon),))
	for i in range(dat.shape[0]):
		interp = intp.RectBivariateSpline(static.x[0,:], static.y[::-1,0], dat[i,::-1,:].T)
		dati[i] = interp.ev(xlon, xlat)

	interp = intp.RectBivariateSpline(static.x[0,:], static.y[::-1,0], ps[::-1,:].T)
	psi = interp.ev(xlon, xlat)
	if not type(p) == type(None):
		interp = intp.RectBivariateSpline(static.x[0,:], static.y[::-1,0], p[::-1,:].T)
		pi = interp.ev(xlon, xlat)
	
	# 2a. Plot the actual map
	xx, xy = m(xlon, xlat)
	if not datmap == None: 
		datmap = __map_prepare_dat(datmap, mask, static, kwargs_map)
		__contourf_dat(m, x, y, datmap, kwargs_map)
	#m.plot(xx, xy, 'w-', linewidth=4)
	m.plot(xx, xy, 'k-', linewidth=4)

	logp = kwargs.get('logp', False)
	if logp:
		z = map(lambda p: np.log(float(p)), static.z)
		kwargs['ticks'] = z
		kwargs['ticklabels'] = static.z
	else:
		z = static.z

	# 2b. Plot the actual cross section
	if 'sect_axes' in kwargs:
		plt.sca(kwargs['sect_axes'])
	else:
		plt.subplot(211)
	
	if kwargs.get('hook'):
		dati = kwargs.pop('hook')(dati)
	xxy = np.array(xxy)
	__contourf_dat(plt, xxy/1e3, z, dati, kwargs)
	if not type(p) == type(None):
		if logp:
			plt.plot(xxy/1e3, np.log(pi), 'g-', linewidth=2)
		else:
			plt.plot(xxy/1e3, pi, 'g-', linewidth=2)
	#plt.gca().invert_yaxis()
	if logp:
		plt.fill_between(xxy/1e3, np.log(psi), np.log(1100.0), color='k')
	else:
		plt.fill_between(xxy/1e3, psi, 1100.0, color='k')
	plt.ylim(float(z[-1]), float(z[0]))
	plt.ylabel('Pressure [hPa]')
	plt.xlabel('Distance along section [km]')

	# 3. Finish off
	# TODO: How to meaningfully replace x,y? Currently the mark feature is broken for sections
	__decorate(m, static, None, slice(None), kwargs)
	__output(kwargs)

	return


def map_oro_dat(dat, static, **kwargs):
	# 1. Prepare
	kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)

	dat = __map_prepare_dat(dat, mask, static, kwargs)
	
	m, x, y, lon, lat = __map_setup(mask, static, kwargs)
	
	# 2. Plot the actual data
	__contourf_dat(m, x, y, dat, kwargs)

	# 3. Finish off
	__decorate(m, x, y, mask, kwargs)
	__output(kwargs)

	return


def map_oro_deform(defabs, defang, static, **kwargs):
	# 1. Prepare
	kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)

	defabs = __map_prepare_dat(defabs, mask, static, kwargs)
	defang = __map_prepare_dat(defang, mask, static, c.contour.defang)
	defdex = np.cos(defang[:,:]) *defabs
	defdey = np.sin(defang[:,:]) *defabs

	# Masking everything below the lowest contour
	if not type(kwargs.get('scale', None)) == type(None):
		defdex[defabs < kwargs.get('scale')[0]] = np.nan
		defdey[defabs < kwargs.get('scale')[0]] = np.nan

	m, x, y, lon, lat = __map_setup(mask, static, kwargs)

	Nvecx = 27
	Nvecy = 18

	#2. Plot the actual deformation
	if hasattr(m, 'transform_vector'):
		ut,vt,xt,yt = m.transform_vector(defdex[::-1,:],defdey[::-1,:],lon[0,:],lat[::-1,0], Nvecx, Nvecy, returnxy=True)
		qscale = 420
	else:
		skipx = defdex.shape[1]/Nvecx
		skipy = defdey.shape[0]/Nvecy
		slc = (slice(skipy/2,None,skipy), slice(skipx/2,None,skipx))
		ut,vt,xt,yt = defdex[slc],defdey[slc],x[slc],y[slc]
		qscale=72
	m.quiver(xt, yt, ut, vt, zorder=4, scale=qscale, alpha=0.85)
	m.quiver(xt, yt, -ut, -vt, zorder=4, scale=qscale, alpha=0.85)
	
	__contourf_dat(m, x, y, defabs, kwargs)
	
	# 3. Finish off
	__decorate(m, x, y, mask, kwargs)
	__output(kwargs)	

	return


def map_oro_barb(u, v, static, dat=None, **kwargs):
	# 1. Prepare
	kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)

	u = __map_prepare_dat(u, mask, static, c.contour.u)
	v = __map_prepare_dat(v, mask, static, c.contour.v)
	if type(dat) == np.ndarray: 
		dat = __map_prepare_dat(dat, mask, static, kwargs)
	
	m, x, y, lon, lat = __map_setup(mask, static, kwargs)
	
	# 2. Plot the actual data
	if type(dat) == np.ndarray: __contourf_dat(m, x, y, dat, kwargs)
	
	# TODO: Generalise, and make the arrow density customisable
	try:
		ut,vt, xt,yt = m.transform_vector(u[::-1,:],v[::-1,:],lon[0,:],lat[::-1,0], 30, 20, returnxy=True)
	except ValueError:
		interval = kwargs.pop('vector_space_interval', 15)
		slc = (slice(interval/2,None,interval), slice(interval/2,None,interval))
		ut,vt, xt,yt = m.rotate_vector(u[slc], v[slc], lon[slc], lat[slc], returnxy=True)
	
	if not kwargs.pop('quiver', False):
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5, zorder=3)
	else:
		m.quiver(xt, yt, ut, vt, zorder=3, scale=kwargs.pop('quiver_length', None), scale_units='width')
	
	# 3. Finish off
	__decorate(m, x, y, mask, kwargs)
	__output(kwargs)

	return


def map_oro_fronts(fronts, froff, static, dat=None, **kwargs):
	# 1. Prepare
	kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)
	
	cfrs = __unflatten_fronts_t(fronts[0], froff[0], minlength=5)
	wfrs = __unflatten_fronts_t(fronts[1], froff[1], minlength=5)
	sfrs = __unflatten_fronts_t(fronts[2], froff[2], minlength=5)

	if not dat == None: 
		dat = __map_prepare_dat(dat, mask, static, kwargs)
	
	m, x, y, lon, lat = __map_setup(mask, static, kwargs)
	
	# 2. Plot the actual data
	if not dat == None: __contourf_dat(m, x, y, dat, kwargs)

	# TODO: Remove conversion from gridpoint indexes to lon/lat once fixed
	for cfr in cfrs:
		lonfr = -180 + cfr[:,0]*0.5
		latfr = 90.0 - cfr[:,1]*0.5
		xfr, yfr = m(lonfr, latfr)
		m.plot(xfr, yfr, 'b-', linewidth=2)
	for wfr in wfrs:
		lonfr = -180 + wfr[:,0]*0.5
		latfr = 90.0 - wfr[:,1]*0.5
		xfr, yfr = m(lonfr, latfr)
		m.plot(xfr, yfr, 'r-', linewidth=2)
	for sfr in sfrs:
		lonfr = -180 + sfr[:,0]*0.5
		latfr = 90.0 - sfr[:,1]*0.5
		xfr, yfr = m(lonfr, latfr)
		m.plot(xfr, yfr, 'm-', linewidth=2)

	# 3. Finish off
	__decorate(m, x, y, mask, kwargs)
	__output(kwargs)

	return

def map_oro_fronts_nc(cfrs, wfrs, sfrs, static, dat=None, **kwargs):
	# 1. Prepare
	kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)
	
	if not dat == None: 
		dat = __map_prepare_dat(dat, mask, static, kwargs)
	
	m, x, y = __map_setup(mask, static, kwargs)
	
	# 2. Plot the actual data
	if not dat == None: __contourf_dat(m, x, y, dat, kwargs)

	for cfr in cfrs:
		lenfr = sum(cfr[:,0] > -200)
		xfr, yfr = m(cfr[:lenfr,1], cfr[:lenfr,0])
		m.plot(xfr, yfr, 'b-', linewidth=2)
	for wfr in wfrs:
		lenfr = sum(wfr[:,0] > -200)
		xfr, yfr = m(wfr[:lenfr,1], wfr[:lenfr,0])
		m.plot(xfr, yfr, 'r-', linewidth=2)
	for sfr in sfrs:
		lenfr = sum(sfr[:,0] > -200)
		xfr, yfr = m(sfr[:lenfr,1], sfr[:lenfr,0])
		m.plot(xfr, yfr, 'm-', linewidth=2)

	# 3. Finish off
	__decorate(m, x, y, mask, kwargs)
	__output(kwargs)

	return


# Helper functions
def __prepare_config(kwargs):
	q = kwargs.pop('q', None)
	if q:
		kwargs = c.contourf.merge(q, **kwargs)
	else:
		kwargs = c.contourf.default.merge(**kwargs)
	
	return kwargs

def __line_prepare_config(kwargs):
	q = kwargs.pop('q', None)
	if q:
		kwargs = c.contour.merge(q, **kwargs)
	else:
		kwargs = c.contour.default.merge(**kwargs)
	
	return kwargs

def __map_create_mask(static, kwargs):
	# Potential override by kwarg
	mask = kwargs.pop('mask', None)
	if type(mask) == np.ndarray:
		return concat1(mask)

	plev = kwargs.pop('plev', None)
	datZ = kwargs.pop('Zdata', None)
	
	if plev and not type(datZ) == np.ndarray:
		f,datZ = metopen(c.file_mstat % {'plev': plev, 'q': 'Z'}, 'mean', cut=c.std_slice[1:], no_static=True)
		if f: f.close()
	if type(datZ) == np.ndarray:
		datZ = concat1(datZ)
		mask = datZ[:,:] < concat1(static.oro[:,:])
	else:
		mask = slice(1,0)
	
	return mask

def __map_prepare_dat(dat, mask, static, kwargs):
	if kwargs.get('hook'):
		dat = kwargs.pop('hook')(dat)
	
	if static.cyclic_ew:
		dat = concat1(dat)
	
	dat[mask] = np.nan

	return dat

def __map_setup(mask, static, kwargs):
	if static.cyclic_ew:
		lon, lat = concat1lonlat(static.x, static.y)
		concat = concat1
	else:
		lon, lat = static.x, static.y
		concat = lambda x: x
	
	m = kwargs.pop('m')
	if not m:
		m = plt
		x, y = lon, lat
	else:
		m = m()
		x, y = m(lon,lat)
	
		lines = m.drawcoastlines(color=kwargs.pop('coastcolor'))
		alpha = kwargs.pop('coast_alpha', None)
		if alpha:
			lines.set_alpha(alpha)

		gridcolor = kwargs.pop('gridcolor')
		if gridcolor:
			meridians = kwargs.pop('meridians', range(0,360,30))
			parallels = kwargs.pop('parallels', range(-75,76,15))
			latmax = kwargs.pop('grid_latmax', parallels[-1])
			alpha = kwargs.pop('grid_alpha', None)
			linestyle = kwargs.pop('grid_linestyle', None)
			dashes = kwargs.pop('grid_dashes', [10, 15])
			
			items = m.drawmeridians(meridians, color=gridcolor, latmax=latmax, dashes=dashes)
			for lines, text in items.values():
				for line in lines:
					if alpha:
						line.set_alpha(alpha)
					if linestyle:
						line.set_linestyle(linestyle)

			items = m.drawparallels(parallels, color=gridcolor, latmax=latmax, dashes=dashes)
			for lines, text in items.values():
				for line in lines:
					if alpha:
						line.set_alpha(alpha)
					if linestyle:
						line.set_linestyle(linestyle)
	
	# Before that the latlon-keyword had no effect :-(
	if basemap_version >= '1.0.7':
		x, y = lon, lat
	
	orocolor = kwargs.pop('orocolor')
	if orocolor:
		m.contour(x, y, concat(static.oro), kwargs.pop('oroscale'), latlon=True, colors=orocolor, 
				alpha=kwargs.pop('oroalpha'), zorder=2)
	if type(mask) == np.ndarray:
		m.contourf(x, y, mask, latlon=True, colors=kwargs.pop('maskcolor'))
	
	return m, x, y, lon, lat

def __contourf_dat(m, x, y, dat, kwargs):
	hatch = kwargs.pop('hatches')
	scale = kwargs.pop('scale')
	if scale == 'auto':
		scale = autoscale(dat, **kwargs)
	cs = m.contourf(x, y, dat, scale, latlon=True, zorder=1, **kwargs)

	# Maximise contrast, by making sure that the last colors of the colorbar 
	# actually are identical to the first/last color in the colormap
	if not type(scale) == int:
		extend = kwargs.get('extend')
		if extend == 'both':
			cs.set_clim(scale[0], scale[-1])
		elif extend == 'max':
			cs.set_clim((scale[0]+scale[1])/2.0, scale[-1])
		elif extend == 'min':
			cs.set_clim(scale[0], (scale[-2]+scale[-1])/2.0)
		else:
			cs.set_clim((scale[0]+scale[1])/2.0, (scale[-2]+scale[-1])/2.0)
	
	return

def __decorate(m, x, y, mask, kwargs):
	if not kwargs.pop('disable_cb'):
		orient = kwargs.pop('cb_orientation', 'vertical')
		spacing = kwargs.pop('cb_tickspacing', 'proportional')
		shrink = kwargs.pop('cb_shrink', 0.8)
		if orient == 'vertical':
			pad = 0.02
			frac = 0.08
		else:
			pad = 0.03
			frac = 0.12
			#pad = 0.02
			#frac = 0.08
		cb = plt.colorbar(ticks=kwargs.pop('ticks'), orientation=orient, 
				shrink=shrink, pad=pad, fraction=frac, spacing=spacing)
		if kwargs.get('ticklabels'): 
			if not orient == 'vertical':
				cb.ax.set_xticklabels(kwargs.pop('ticklabels'))
			else:
				cb.ax.set_yticklabels(kwargs.pop('ticklabels'))
	
	#legend_labels = kwargs.pop('legend_labels', None)
	#if legend_labels:
	#	plt.legend([], legend_labels)
	
	for overlay in kwargs.pop('overlays'):
		overlay(m,x,y, zorder=3, mask=mask)
	
	if kwargs.get('mark'):
		yidx, xidx = kwargs.pop('mark')
		m.scatter(x[yidx,xidx], y[yidx,xidx], 484, latlon=True, marker='o', facecolors=(1,1,1,0), 
				edgecolors='k', linewidths=3, zorder=3)
	
	if kwargs.get('title'):
		plt.title(kwargs.pop('title'))
	
	return

def __output(kwargs):
	if kwargs.get('save'):
		plt.savefig(kwargs.pop('save')) #, format='png')
	if kwargs.pop('show'):
		plt.show()
	
	return


# Overlays
def sect_overlay_dat(dat, sect, **kwargs):  
	kwargs = __line_prepare_config(kwargs)

	def overlay(m, static, dump, zorder, mask=None):
		xlon, xlat, xxy = sect_gen_points(sect, m, 25000.0)
		xxy = np.array(xxy)

		dati = np.empty((dat.shape[0], len(xlon),))
		for i in range(dat.shape[0]):
			interp = intp.RectBivariateSpline(static.x[0,:], static.y[::-1,0], dat[i,::-1,:].T)
			dati[i] = interp.ev(xlon, xlat)
	
		if type(mask) == np.ndarray:
			dati[mask] = np.nan
		scale = kwargs.pop('scale')
		if scale == 'auto':
			scale = autoscale(dat, **kwargs)
		if kwargs.get('logp', False):
			z = map(lambda x: np.log(float(z)), static.z)
		else:
			z = static.z
		cs =  plt.contour(xxy/1e3, static.z, dati, scale, **kwargs)

		labels = kwargs.pop('contour_labels')
		if labels:
			plt.clabel(cs, kwargs.pop('contour_labels_scale', scale),
					fontsize=kwargs.pop('contour_labels_fontsize', 12), 
					fmt=kwargs.pop('contour_labels_format', '%1.1f'),
					inline=kwargs.pop('contour_labels_inline', True), 
					inline_spacing=kwargs.pop('contour_labels_inline_spacing', 2)
			)

		return

	return overlay


def map_overlay_dat(dat, static, **kwargs):  
	kwargs = __line_prepare_config(kwargs)
	
	if static.cyclic_ew:
		dat = concat1(dat)
	if kwargs.get('hook'):
		dat = kwargs.pop('hook')(dat)

	def overlay(m, x, y, zorder, mask=None):
		if type(mask) == np.ndarray:
			dat[mask] = np.nan
		scale = kwargs.pop('scale')
		if scale == 'auto':
			scale = autoscale(dat, **kwargs)
		cs =  m.contour(x, y, dat, scale, latlon=True, **kwargs)

		labels = kwargs.pop('contour_labels')
		if labels:
			plt.clabel(cs, kwargs.pop('contour_labels_scale', scale),
					fontsize=kwargs.pop('contour_labels_fontsize', 12), 
					fmt=kwargs.pop('contour_labels_format', '%1.1f'),
					inline=kwargs.pop('contour_labels_inline', True), 
					inline_spacing=kwargs.pop('contour_labels_inline_spacing', 2)
			)

		return

	return overlay


def map_overlay_fronts(fronts, froff, static, **kwargs):  
	kwargs = __line_prepare_config(kwargs)

	cfrs = __unflatten_fronts_t(fronts[0], froff[0], minlength=5)
	wfrs = __unflatten_fronts_t(fronts[1], froff[1], minlength=5)
	sfrs = __unflatten_fronts_t(fronts[2], froff[2], minlength=5)

	def overlay(m, x, y, zorder, mask=None):
		# TODO: Remove conversion from gridpoint indexes to lon/lat once fixed
		# TODO: Convert to latlon=True system
		for cfr in cfrs:
			if hasattr(m, '__call__'):
				lonfr = static.x[0,0] + (cfr[:,0]-1)*(static.x[0,1]-static.x[0,0])
				latfr = static.y[0,0] + (cfr[:,1]-1)*(static.y[1,0]-static.y[0,0])
				xfr, yfr = m(lonfr, latfr)
			else:
				xfr, yfr = cfr[:,0]*50000, cfr[:,1]*50000
			m.plot(xfr, yfr, 'b-', linewidth=2)
		for wfr in wfrs:
			if hasattr(m, '__call__'):
				lonfr = static.x[0,0] + (wfr[:,0]-1)*(static.x[0,1]-static.x[0,0])
				latfr = static.y[0,0] + (wfr[:,1]-1)*(static.y[1,0]-static.y[0,0])
				xfr, yfr = m(lonfr, latfr)
			else:
				xfr, yfr = wfr[:,0]*50000, wfr[:,1]*50000
			m.plot(xfr, yfr, 'r-', linewidth=2)
		for sfr in sfrs:
			if hasattr(m, '__call__'):
				lonfr = static.x[0,0] + (sfr[:,0]-1)*(static.x[0,1]-static.x[0,0])
				latfr = static.y[0,0] + (sfr[:,1]-1)*(static.y[1,0]-static.y[0,0])
				xfr, yfr = m(lonfr, latfr)
			else:
				xfr, yfr = sfr[:,0]*50000, sfr[:,1]*50000
			m.plot(xfr, yfr, 'm-', linewidth=2)


		return

	return overlay


def map_overlay_lines(lines, loff, static, **kwargs):  
	kwargs = __line_prepare_config(kwargs)

	lns = __unflatten_fronts_t(lines, loff, minlength=0)

	def overlay(m, x, y, zorder, mask=None):
		# TODO: Remove conversion from gridpoint indexes to lon/lat once fixed
		# TODO: Convert to latlon=True system
		for ln in lns:
			lonfr = static.x[0,0] + (ln[:,0]-1)*(static.x[0,1]-static.x[0,0])
			latfr = static.y[0,0] + (ln[:,1]-1)*(static.y[1,0]-static.y[0,0])
			xfr, yfr = m(lonfr, latfr)
			m.plot(xfr, yfr, kwargs['linecolor'], linewidth=kwargs.get('linewidth', 2), 
					alpha=kwargs.get('alpha', 1))


		return

	return overlay


def map_overlay_dots(xidxs, yidxs, static, **kwargs):  
	kwargs = __line_prepare_config(kwargs)

	def overlay(m, x, y, zorder, mask=None):
		# TODO: Remove conversion from gridpoint indexes to lon/lat once fixed
		# TODO: Convert to latlon=True system
		lonfr = static.x[0,0] + (xidxs - 1)*(static.x[0,1]-static.x[0,0])
		latfr = static.y[0,0] + (yidxs - 1)*(static.y[1,0]-static.y[0,0])
		xfr, yfr = m(lonfr, latfr)
		m.scatter(xfr, yfr, kwargs.get('linewidth',9), marker='.', edgecolors=kwargs['linecolor'])


		return

	return overlay


def map_overlay_mask(dat, **kwargs):  
	kwargs = __line_prepare_config(kwargs)

	dat = concat1(dat)

	def overlay(m, x, y, zorder, mask=None):
		if type(mask) == np.ndarray:
			raise NotImplementedError, 'here be dragons'
		cs =  m.contourf(x, y, dat, latlon=True, **kwargs)

		return

	return overlay


def map_overlay_latlonbox(lon0, lon1, lat0, lat1, vertices=30, **kwargs):
	def overlay(m, x, y, zorder, mask=None):
		# Western boundary
		m.plot(np.ones((vertices,))*lon0, np.linspace(lat0,lat1,vertices), latlon=True, **kwargs)
		# Eastern boundary
		m.plot(np.ones((vertices,))*lon1, np.linspace(lat0,lat1,vertices), latlon=True, **kwargs)

		# Southern boundary
		m.plot(np.linspace(lon0,lon1,vertices), np.ones((vertices,))*lat0, latlon=True, **kwargs)
		# Northern boundary
		m.plot(np.linspace(lon0,lon1,vertices), np.ones((vertices,))*lat1, latlon=True, **kwargs)

		return

	return overlay


def phist_defang(defang, defabs=None):
	pass


# that's it
