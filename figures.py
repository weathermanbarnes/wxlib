#!/usr/bin/python
# -*- encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm as bmcm
from metopen import metopen
from utils import concat1

import static as c


# globally useful
f = metopen('static.npz')
oro = f['oro']
oro = concat1(oro)
s   = oro.shape
lat = np.tile(f['lat'], (s[1],1)).T
lon = np.tile(f['lon'], (s[0],1))
lon = concat1(lon)
f.close()

orolevs = range(-19000,51000,2000)

# 32 years of data at a glance
def ypline_mean(q='defabs', yidx=57, plev=700, summarize=False, std=False, quiet=False):
	if not quiet:
		print 'Lat %f' % lat[yidx,0]
	means = {}
	for y in c.years:
		if not quiet:
			print y
		f = metopen(c.file_stat % (y, plev, q))
		if std:
			means[y] = f['stddev'][yidx,:]
		else:
			means[y] = f['mean'][yidx,:]
		f.close()
	
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
	meanZ = np.ones(s)*9999999
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

		f   = metopen(c.file_stat % (y, plev, q))
		#fZ  = metopen(c.file_stat % (y, plev, 'Z'))
		if std:	
			# TODO: Correct calculation of mean stddev!
			mean += concat1(f['stddev'])
		else:
			mean += concat1(f['mean'])
		#meanZ += concat1(fZ['mean'])

		f.close()
		#fZ.close()
	
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
		fileu = 'mean_u.npz'
		filev = 'mean_v.npz'
		fileZ = 'mean_Z.npz'
	else: 
		fileu = '%04d_umean.npz' % year
		filev = '%04d_vmean.npz' % year
		fileZ = '%04d_Zmean.npz' % year

	f = metopen(fileu)
	meanu = f['u%d' % plev]
	f.close()

	f = metopen(filev)
	meanv = f['v%d' % plev]
	f.close()
	
	f = metopen(fileZ)
	meanZ = f['Z%d' % plev]
	f.close()

	f = pygrib.open('static.grib')
	oro = f[3]
	lats,lons = oro.latlons()
	f.close()

	lats = np.concatenate((lats, np.reshape(lats[:,0], (41,1)) ), axis=1)
	lons = np.concatenate((lons, np.ones((41,1))*179.999 ), axis=1)
	oro = np.concatenate((oro.values, np.reshape(oro.values[:,0], (41,1)) ), axis=1)
	meanu= np.concatenate((meanu, np.reshape(meanu[:,0], (41,1)) ), axis=1)
	meanv= np.concatenate((meanv, np.reshape(meanv[:,0], (41,1)) ), axis=1)
	meanZ= np.concatenate((meanZ, np.reshape(meanZ[:,0], (41,1)) ), axis=1)
	
	meanu = 3*meanu
	meanv = 3*meanv
	print meanu.min(), meanu.max(), meanu.mean()
	print meanv.min(), meanv.max(), meanv.mean()
	meanu[meanZ < oro] = np.ma.masked
	meanv[meanZ < oro] = np.ma.masked

	m = Basemap(projection='npstere',boundinglat=30,lon_0=-50,resolution='l')
	x,y = m(lons,lats)
	ut,vt,xt,yt = m.transform_vector(meanu,meanv,lons[0,:],lats[:,0],120,120,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth)
	if not quiver:
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5)
	else:
		m.quiver(xt, yt, ut, vt)
	m.drawparallels(range(30,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# Contour map of averaged wind on top of a oro map.
def map_mean_deformv(year=None, plev=700, quiver=True):
	if not year:
		fileFF = 'mean_rotdef.npz'
		fileDD = 'mean_defang.npz'
		fileZ = 'mean_Z.npz'
		nameext = 'm'
	else: 
		fileFF = '%04d_rotdef.npz' % year
		fileDD = '%04d_defang.npz' % year
		fileZ = '%04d_Zmean.npz' % year
		nameext = '%02d' % (year % 100)

	f = metopen(fileFF)
	meanFF = f['d%s_%d' % (nameext, plev)]
	f.close()

	f = metopen(fileDD)
	meanDD = f['da%s_%04d' % (nameext, plev)]
	f.close()
	
	#if year:
	#	meanDD = meanDD.mean(axis=0)
	#	meanFF = meanFF.mean(axis=0)
	
	meanu = np.cos(meanDD[:,:]) *meanFF
	meanv = np.sin(meanDD[:,:]) *meanFF
	
	f = metopen(fileZ)
	meanZ = f['Z%d' % plev]
	f.close()

	f = pygrib.open('static.grib')
	oro = f[3]
	lats,lons = oro.latlons()
	f.close()

	if year:
		meanu = meanu.mean(axis=0)
		meanv = meanv.mean(axis=0)
		meanFF= meanFF.mean(axis=0)
	#print meanZ.shape

	lats  = np.concatenate((lats, np.reshape(lats[:,0], (41,1)) ), axis=1)
	lons  = np.concatenate((lons, np.ones((41,1))*179.999 ), axis=1)
	oro  = np.concatenate((oro.values, np.reshape(oro.values[:,0], (41,1)) ), axis=1)
	meanu = np.concatenate((meanu, np.reshape(meanu[:,0], (41,1)) ), axis=1)
	meanv = np.concatenate((meanv, np.reshape(meanv[:,0], (41,1)) ), axis=1)
	meanFF= np.concatenate((meanFF, np.reshape(meanFF[:,0], (41,1)) ), axis=1) 
	meanZ = np.concatenate((meanZ, np.reshape(meanZ[:,0], (41,1)) ), axis=1)
	
	print meanu.min(), meanu.max(), meanu.mean()
	print meanv.min(), meanv.max(), meanv.mean()
	meanu [meanZ < oro] = np.ma.masked
	meanv [meanZ < oro] = np.ma.masked
	meanFF[meanZ < oro] = np.ma.masked

	m = Basemap(projection='npstere',boundinglat=30,lon_0=-50,resolution='l')
	x,y = m(lons,lats)
	ut,vt,xt,yt = m.transform_vector(meanu,meanv,lons[0,:],lats[:,0],60,60,returnxy=True)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth, zorder=1)
	m.contour(x[1:-2], y[1:-2], meanFF[1:-2], 25, zorder=2)
	if not quiver:
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5, zorder=3)
	else:
		m.quiver(xt, yt, ut, vt, zorder=3)
		m.quiver(xt, yt,-ut,-vt, zorder=3)
	m.drawparallels(range(30,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def map_date_deform(date, plev=700, norm=False, std=False):
	tidx = (date.timetuple().tm_yday-1)*4 + int(date.hour/6)
	
	if norm:
		f = metopen('%04d_deforn.npz' % date.year)
	else:
		f = metopen('%04d_deform.npz' % date.year)
	
	d = f['d%02d_%d' % (date.year % 100, plev)][tidx,:,:]
	f.close()

	f = pygrib.open('%04d_Z.grib' % date.year)
	Z = f[13*tidx+list(levels).index(plev)+1]
	print Z.year, Z.month, Z.day, Z.hour, Z.level
	Z = Z.values
	f.close()
	
	f = pygrib.open('static.grib')
	oro = f[3]
	lats,lons = oro.latlons()
	f.close()

	lats = np.concatenate((lats, np.reshape(lats[:,0], (41,1)) ), axis=1)
	lons = np.concatenate((lons, np.ones((41,1))*179.999 ), axis=1)
	oro = np.concatenate((oro.values, np.reshape(oro.values[:,0], (41,1)) ), axis=1)
	d    = np.concatenate((d, np.reshape(d[:,0], (39,1)) ), axis=1)
	Z    = np.concatenate((Z, np.reshape(Z[:,0], (41,1)) ), axis=1)

	print oro.shape, levels.shape
	d[Z < oro] = np.ma.masked

	m = Basemap(projection='npstere',boundinglat=30,lon_0=-50,resolution='l')
	x,y = m(lons,lats)
	m.drawcoastlines()
	m.contourf(x, y, oro, orolevs, cmap=plt.cm.gist_earth)
	m.contour(x[1:-1,:], y[1:-1,:], d, 20)
	m.drawparallels(range(30,80,5))
	m.drawmeridians(range(0,360,30))
	plt.colorbar()
	plt.show()

	return


# vertical profiles of 32years mean deformation
def ysect_mean_deform(year=None, yidx=57, rot=False, norm=False, filename=None, std=False, quiet=False):
	dyidx = 0
	if not year:
		if filename:
			f = metopen(filename+'.npz')
		if norm:
			f = metopen('mean_deforn.npz')
		elif rot:
			f = metopen('mean_rotdef.npz')
			dyidx = 1
		else:
			f = metopen('mean_deform.npz')
		ln = 3
	else:
		if filename:
			f = metopen(filename+'.npz')
		elif norm:
			f = metopen('%04d_deforn.npz' % year)
		elif rot:
			f = metopen('%04d_rotdef.npz' % year)
			dyidx = 1
		else:
			f = metopen('%04d_defor2.npz' % year)
		ln = 4

	dm = np.zeros((13,240))
	i = 0
	for field in sorted(map(lambda x: int(x[ln:]), f.files)):
		if not year:
			dm[i] = f['dm_%d' % field][yidx+dyidx,:]
		else :
			dm[i] = f['d%02d_%d' % (year % 100, field)][:,yidx+dyidx,:].mean(axis=0)
		i += 1
	
	fs = pygrib.open('static.grib')
	oro = fs[3]
	fs.close()
	#oro = 1000.0-oro.values[yidx,:]/80.0	
	oro = 1000.0*np.exp(oro.values[yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	oro[oro > 1000.0] = 1000.0

	print oro.shape, levels.shape
	dm[oro[np.newaxis,:]-10 < levels[:,np.newaxis]] = np.ma.masked

	plt.contour(lon[0,:], levels, dm, 20)
	#plt.fill_between(lon[0,:], 1000.0-oro.values[yidx,:]/80.0, np.ones((240,))*1000.0, 'k')
	plt.fill(lon[0,:], oro, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


# vertical profiles of 32years mean deformation
def ysect_meanQ(quantity, yidx=57, norm=False, std=False, quiet=False):
	f = metopen('mean_%s.npz' % quantity)

	Qm = np.zeros((13,240))
	i = 0
	for field in sorted(map(lambda x: int(x[len(quantity):]), f.files)):
		Qm[i] = f['%s%d' % (quantity, field)][yidx,:]
		i += 1
	
	fs = pygrib.open('static.grib')
	oro = fs[3]
	fs.close()
	#oro = 1000.0-oro.values[yidx,:]/80.0	
	oro = 1000.0*np.exp(oro.values[yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	oro[oro > 1000.0] = 1000.0

	print oro.shape, levels.shape
	Qm[oro[np.newaxis,:]-10 < levels[:,np.newaxis]] = np.ma.masked

	plt.contour(lon[0,:], levels, Qm, 15)
	#plt.fill_between(lon[0,:], 1000.0-oro.values[yidx,:]/80.0, np.ones((240,))*1000.0, 'k')
	plt.fill(lon[0,:], oro, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def ysect_date_deform(date, yidx=57, norm=False, std=False):
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


# same as ysect_mean_deform but without any averaging; deformation sections for one point in time
def ysect_dateQ(quantity, date, yidx=57, norm=False, std=False):
	tidx = (date.timetuple().tm_yday-1)*4 + int(date.hour/6)
	
	f = pygrib.open('%04d_%s.grib' % (date.year, quantity))

	Q = np.zeros((13,240))
	for i in range(13):
		print levels[i], f[13*tidx+i+1].level
		Q[i] = f[13*tidx+i+1].values[yidx,:]
	
	fs = pygrib.open('static.grib')
	oro = fs[3]
	fs.close()
	oro = 1000.0*np.exp(oro.values[yidx,:]/(-270.0*287.0))	# tentative conversion from Phi [m^2/s^2] to p [hPa] p0 = 1000, <T> = 270K
	oro[oro > 1000.0] = 1000.0

	print oro.shape, levels.shape
	Q[oro[np.newaxis,:]-10 < levels[:,np.newaxis]] = np.ma.masked

	plt.contour(lon[0,:], levels, Q, 15)
	plt.fill(lon[0,:], oro, 'k')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


# Hovmöller diagram
def ypline_hov(year, plev=700, yidx=57, norm=False):
	if norm:
		f = metopen('%04d_deforn.npz' % year)
	else:
		f = metopen('%04d_deform.npz' % year)

	dat = f['d%02d_%d' % (year % 100, plev)]

	dat = np.concatenate((dat[:,:,-30:], dat, dat[:,:,:30]), axis=2)
	exlon = map(lambda x: x*1.5 - 225.0, range(dat.shape[2]))
	tidxs = map(lambda x: x*0.25,range(dat.shape[0]))

	plt.contourf(exlon, tidxs, dat[:,yidx,:], 20)
	plt.plot([-180, -180], [tidxs[0], tidxs[-1]], 'k--')
	plt.plot([ 180,  180], [tidxs[0], tidxs[-1]], 'k--')
	plt.ylim(plt.ylim()[::-1])		# reverse y-axis
	plt.colorbar()
	plt.show()

	return


def correl_mean_deform_v():
	fv = metopen('mean_v.npz')
	fd = metopen('mean_deform.npz')
	
	v = []
	for f in fv.files:
		v.append(fv[f][1:-1,:])
	d = []
	for f in fd.files:
		d.append(fd[f])

	plt.scatter(d, v, s=1, marker='+')
	plt.show()

	return


def xsect_mean_deform(year=None, norm=False):
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



