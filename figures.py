#!/usr/bin/python
# -*- encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm as bmcm
from metopen import metopen
from utils import concat1

import static as c


# globally useful
f, oro = metopen('static', 'oro')
oro = concat1(oro)
s   = oro.shape
lat = np.tile(f['lat'], (s[1],1)).T
lon = np.tile(f['lon'], (s[0],1))
lon = concat1(lon)
lon[:,-1] += 360
f.close()

orolevs = range(-19000,51000,2000)

# 32 years of data at a glance
def ypline_mean_Q(q='defabs', yidx=57, plev=700, summarize=False, std=False, quiet=False):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	means = {}
	for y in c.years:
		if not quiet:
			print y
		if std:
			f, dat = metopen(c.file_stat % (y, plev, q), 'stddev')
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
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# contour map of 32 year mean deformation
def map_mean_Q(q='defabs', year=None, plev=700, std=False, quiet=False):
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
		if not quiet:
			print y

		if std:	
			# TODO: Correct calculation of mean stddev!
			f, dat = metopen(c.file_stat % (y, plev, q), 'stddev')
		else:
			f, dat = metopen(c.file_stat % (y, plev, q), 'mean')
		mean = concat1(dat)

		fZ, daZ = metopen(c.file_stat % (y, plev, 'Z'), 'mean')
		meanZ += concat1(daZ)
		
		del dat, daZ
		f.close()
		fZ.close()
	
	mean  /= len(years)
	meanZ /= len(years)
	
	# Mask out surface parts below orography
	mean[meanZ < oro] = np.ma.masked
	mean[0:5,:] = np.ma.masked

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	m.drawcoastlines()
	x,y = m(lon,lat)
	m.contourf(x,y,oro,orolevs, cmap=plt.cm.gist_earth)
	m.contour(x, y, mean, 20)
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Contour map of averaged wind on top of a oro map.
def map_mean_barb(year=None, plev=700, quiver=False):
	if not year:
		fu, meanu = metopen(c.file_mstat % (plev, 'u'), 'mean')
		fv, meanv = metopen(c.file_mstat % (plev, 'v'), 'mean')
		fZ, meanZ = metopen(c.file_mstat % (plev, 'Z'), 'mean')
	else: 
		fu, meanu = metopen(c.file_stat % (year, plev, 'u'), 'mean')
		fv, meanv = metopen(c.file_stat % (year, plev, 'v'), 'mean')
		fZ, meanZ = metopen(c.file_stat % (year, plev, 'Z'), 'mean')

	meanu = concat1(meanu[::-1,:]*5)
	meanv = concat1(meanv[::-1,:]*5)
	meanZ = concat1(meanZ[::-1,:])
	
	print meanu.min(), meanu.max(), meanu.mean()
	print meanv.min(), meanv.max(), meanv.mean()
	meanu[meanZ < oro] = np.ma.masked
	meanv[meanZ < oro] = np.ma.masked

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


# Contour map of averaged deformation vector on top of a oro map.
def map_mean_deform(year=None, plev=700):
	if not year:
		fabs, meanabs = metopen(c.file_mstat % (plev, 'defabs'), 'mean')
		fang, meanang = metopen(c.file_mstat % (plev, 'defang'), 'mean')
		fZ, meanZ     = metopen(c.file_mstat % (plev, 'Z'), 'mean')
	else: 
		fabs, meanabs = metopen(c.file_stat % (year, plev, 'defabs'), 'mean')
		fang, meanang = metopen(c.file_stat % (year, plev, 'defang'), 'mean')
		fZ, meanZ     = metopen(c.file_stat % (year, plev, 'Z'), 'mean')
	
	meandex = np.cos(meanang[:,:]) *meanabs
	meandey = np.sin(meanang[:,:]) *meanabs
	
	meandex = concat1(meandex[::-1,:])
	meandey = concat1(meandey[::-1,:]) 
	meanabs = concat1(meanabs)
	meanZ   = concat1(meanZ)
	
	print meandex.min(), meandex.max(), meandex.mean()
	print meandey.min(), meandey.max(), meandey.mean()
	meanabs[meanZ < oro] = np.ma.masked
	meanZ = meanZ[::-1,:]
	meandex[meanZ < oro] = np.ma.masked
	meandey[meanZ < oro] = np.ma.masked

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	ut,vt,xt,yt = m.transform_vector(meandex,meandey,lon[0,:],lat[::-1,0],60,60,returnxy=True)
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


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def map_date_Q(date, q='defabs', plev=700):
	tidx = (date.timetuple().tm_yday-1)*4 + int(date.hour/6)
	
	f, dat = metopen(c.file_std % (date.year, plev, q), c.q[q])
	dat = dat[tidx,:,:]
	f.close()

	f, daZ = metopen(c.file_std % (date.year, plev, 'Z'), c.q['Z'])
	daZ = daZ[tidx,:,:]
	f.close()

	dat = concat1(dat)
	daZ = concat1(daZ)

	dat[daZ < oro] = np.ma.masked

	m = Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
	x,y = m(lon,lat)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth)
	m.contour(x, y, dat, 25)
	m.drawparallels(range(15,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# vertical profiles of 32years mean deformation
def ysect_mean_Q(q='defabs', yidx=57, std=False, quiet=False):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	i = 0
	qm = np.zeros((len(c.plevs),s[1]))
	for plev in c.plevs:
		if std:
			f, dat = metopen(c.file_mstat % (plev, q), 'stddev')
		else:
			f, dat = metopen(c.file_mstat % (plev, q), 'mean')

		qm[i] = concat1(dat[yidx,:])
		i += 1
	
	oro = 1000.0*np.exp(globals()['oro'][yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	oro[oro > 1000.0] = 1000.0
	
	qm[(oro-10) < c.plevs[:,np.newaxis]] = np.ma.masked

	plt.contour(lon[0,:], c.plevs, qm, 15)
	#plt.fill_between(lon[0,:], 1000.0-oro.values[yidx,:]/80.0, np.ones((240,))*1000.0, 'k')
	plt.fill(lon[0,:], oro, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def ysect_date_Q(date, yidx=57, norm=False, std=False):
	raise NotImplementedError, "This would be _very_ slow. Make it fast and remove this error!"
	
	tidx = (date.timetuple().tm_yday-1)*4 + int(date.hour/6)
	
	if norm:
		f = metopen('%04d_deforn.npz' % date.year)
	else:
		f = metopen('%04d_deform.npz' % date.year)
	d = np.zeros((13,240))
	i = 0
	for field in sorted(map(lambda x: int(x[4:]), f.files)):
		d[i] = f['d%02d_%d' % (date.year % 100, field)][tidx,yidx,:]
		i += 1
	
	fs = pygrib.open('static.grib')
	oro = fs[3]
	fs.close()
	oro = 1000.0*np.exp(oro.values[yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	oro[oro > 1000.0] = 1000.0

	print oro.shape, levels.shape
	d[oro[np.newaxis,:]-10 < levels[:,np.newaxis]] = np.ma.masked

	plt.contour(lon[0,:], levels, d, 15)
	plt.fill(lon[0,:], oro, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


# Hovmöller diagram
def ypline_hov_Q(year, q='defabs', plev=700, yidx=57, norm=False, quiet=False):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	f, dat = metopen(c.file_std % (year, plev, q), c.q[q])
	dat = dat[:,yidx,:]
	f.close()

	dat = np.concatenate((dat[:,-60:], dat, dat[:,:60]), axis=1)
	exlon = map(lambda x: x*0.5 - 210.0, range(dat.shape[1]))
	tidxs = map(lambda x: x*0.25,range(dat.shape[0]))

	plt.contourf(exlon, tidxs, dat[:,:], 20)
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


def xsect_mean_Q(q='defabs', year=None, norm=False):
	if not year:
		if norm:
			f = metopen('mean_deforn.npz')
		else:
			f = metopen('mean_deform.npz')
		ln = 3
	else:
		if norm:
			f = metopen('%04d_deforn.npz' % year)
		else:
			f = metopen('%04d_deform.npz' % year)
		ln = 4

	dm = np.zeros((13,39))
	i = 0
	for field in sorted(map(lambda x: int(x[ln:]), f.files)):
		if not year:
			dm[i] = f['dm_%d' % field].mean(axis=1)
		else:
			dm[i] = f['d%02d_%d' % (year % 100, field)].mean(axis=0).mean(axis=1)
		i += 1
	
	plt.contour(lat[1:-1,0], levels, dm, 15)
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return



