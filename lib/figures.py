#!/usr/bin/env python
# -*- encoding: utf-8

''' A collection of general plotting functions

These functions have two main aims:

 1. Work around some of the peculiarities and sharp edges of matplotlib.
    This includes making some often-used plot customisation items more easily available 
    by keyword argument, and the automatic concatenation of the 180°E-Meridian to represent
    the periodicity of the data in global contour(f) plots.
 2. Take into account a large repository of default plot configuration, such that
    for example temperatures are by default plotted in a consistent an meaningful way.
'''


import copy

import numpy as np
import scipy.interpolate as intp

import matplotlib.pyplot as plt
import mpl_toolkits.basemap
basemap_version = mpl_toolkits.basemap.__version__

from settings import conf

from metio import metopen
from utils import concat1, concat1lonlat, __unflatten_fronts_t, sect_gen_points
from autoscale import autoscale


def section_p(dat, ps, sect, static, datmap=None, p=None, **kwargs):
	''' Plot a vertical cross section with pressure or log-pressure as the vertical axis

	The plot is highly customisable. It features an inset map showing the location 
	of the cross section on a map. The cross section follows a straight line in the 
	given map projection. Using different map projections one can hence construct cross 
	sections along the following special lines:

	 * Great circle: Use the gnomic projection
	 * Latitude circle: Use a cylindrical projection or plain latlon coordinates
	 * Longitude circle: Use the gnomic projection or plain latlon coordinates
	 * Rhumb lines (loxodromes or lines of constant course): Use the Mercator 
	   projection
	
	For other projections, straight lines generally do not have special properties.

	The data is interpolated along this straight line using equidistant 
	interpolation points. By default the interpolation point distance is about
	25 km, appropriate for cross sections covering larger regions.

	Parameters
	----------
	dat : np.ndarray with the dimensions (z,y,x)
		Data array from which to take the cross section
	ps : np.ndarray with the dimensions (y,x)
		Surface pressure for masking areas below ground
	sect : list of 2-tuples
		List of coordinates defining the cross section
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	datmap : np.ndarray with the dimensions (y,x)
		*Optional*: Data to be displayed on the inset map. 
	p : np.ndarray with the dimensions (y,x)
		*Optional*: Mark the height of a given surface on the cross section
		by supplying its height in pressure coordinates.
	
	Keyword arguments
	-----------------
	plot arguments : all contourf
		For a list of valid arguments refer to :ref:`plot configuration`. Keyword 
		arguments are applied both to the inset map and the main cross section.
	'''
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
	# TODO: Make the 25 km configurable
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
	
	# 2a. Plot the inset map
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

# Section using pressure as vertical coordinate are standard, so provide shortcut
section = section_p


def map(dat, static, **kwargs):
	''' Plot data on a map

	The plot is highly customisable and takes dynlib defaults as well as 
	personal settings into account.

	Parameters
	----------
	dat : np.ndarray with dimensions (y,x)
		Data to be plotted
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contourf
		For a list of valid arguments refer to :ref:`plot configuration`.
	'''
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


# TODO: The dilatation axes should become an overlay, and this function should be removed
def map_deform(defabs, defang, static, **kwargs):
	''' Soon obsolete '''
	# 1. Prepare
	kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)

	defabs = __map_prepare_dat(defabs, mask, static, kwargs)
	defang = __map_prepare_dat(defang, mask, static, conf.contour.defang)
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


# TODO: Wind barbs and arrows should become overlays, and this function should be removed
def map_barb(u, v, static, dat=None, **kwargs):
	''' Soon obsolete '''
	# 1. Prepare
	kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)

	u = __map_prepare_dat(u, mask, static, conf.contour.u)
	v = __map_prepare_dat(v, mask, static, conf.contour.v)
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


# TODO: Fronts should become an overlay, and this function should be removed
# TODO: Merge with other fronts overlay
def map_oro_fronts(fronts, froff, static, dat=None, **kwargs):
	''' Obsolete: Use map_overlay_fronts instead '''
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

# TODO: Fronts should become an overlay, and this function should be removed
# TODO: Merge with other fronts overlay
def map_oro_fronts_nc(cfrs, wfrs, sfrs, static, dat=None, **kwargs):
	''' Obsolete: Use map_overlay_fronts instead '''
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
	''' Make sure kwargs contains a complete contourf plot configuration, 
	filling undefined keys from the dynlib configuration.'''
	q = kwargs.pop('q', None)
	if q:
		kwargs = conf.contourf.merge(q, **kwargs)
	else:
		kwargs = conf.contourf.default.merge(**kwargs)
	
	return kwargs

def __line_prepare_config(kwargs):
	''' Make sure kwargs contains a complete contour plot configuration, 
	filling undefined keys from the dynlib configuration.'''
	q = kwargs.pop('q', None)
	if q:
		kwargs = conf.contour.merge(q, **kwargs)
	else:
		kwargs = conf.contour.default.merge(**kwargs)
	
	return kwargs

def __map_create_mask(static, kwargs):
	''' Consider masked areas
	
	 1. Explicitly given 
	 2. Orography using instantaneous z
	 3. Orography using climatological z
	'''
	# Potential override by kwarg
	mask = kwargs.pop('mask', None)
	if type(mask) == np.ndarray:
		return concat1(mask)

	plev = kwargs.pop('plev', None)
	datZ = kwargs.pop('Zdata', None)
	
	if plev and not type(datZ) == np.ndarray:
		f,datZ = metopen(conf.file_mstat % {'plev': plev, 'q': 'Z'}, 'mean', cut=conf.std_slice[1:], no_static=True)
		if f: f.close()
	if type(datZ) == np.ndarray:
		datZ = concat1(datZ)
		mask = datZ[:,:] < concat1(static.oro[:,:])
	else:
		mask = slice(1,0)
	
	return mask

def __map_prepare_dat(dat, mask, static, kwargs):
	''' Prepare the data to be plotted '''
	if kwargs.get('hook'):
		dat = kwargs.pop('hook')(dat)
	
	if static.cyclic_ew:
		dat = concat1(dat)
	
	dat[mask] = np.nan

	return dat

def __map_setup(mask, static, kwargs):
	''' Setup the map projection and draw the lat/lon grid '''
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
	''' Plot the actual data '''
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
	''' Add "decorations": colorbar, legends, overlays and a title'''
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
	''' Save and or show the plot '''
	if kwargs.get('save'):
		plt.savefig(kwargs.pop('save'))
	if kwargs.pop('show'):
		plt.show()
	
	return


###############################################
# Overlays

# TODO: Should this one take the static object as an argument? If only for consistency in the API.
def section_overlay_dat(dat, sect, **kwargs):
	''' Overlay contours onto a section
	
	Parameters
	----------
	dat : np.ndarray with dimensions (z,y,x)
		Data to be plotted
	sect : list of 2-tuples
		List of coordinates defining the cross section
	
	Keyword arguments
	-----------------
	plot arguments : all contour
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
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
	''' Overlay contours onto a map
	
	Parameters
	----------
	dat : np.ndarray with dimensions (z,y,x)
		Data to be plotted
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contour
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
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
	''' Overlay front lines onto a map
	
	Parameters
	----------
	fronts : np.ndarray with dimensions (fronttype,pointindex,infotype)
		Front position array
	fronts : np.ndarray with dimensions (fronttype,frontindex)
		Front offset array, a list of starting point indexes
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contour
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
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
	''' Overlay lines onto a map
	
	Parameters
	----------
	lines : np.ndarray with dimensions (pointindex,infotype)
		Line position array
	loff : np.ndarray with dimension (lineindex)
		Line offset array, a list of starting point indexes
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contour
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
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
	''' Overlay dots onto a map
	
	Parameters
	----------
	xidxs, yidxs : np.ndarrays with dimension (pointindex)
		Dot position lists
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contour
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
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


def map_overlay_shading(dat, static, **kwargs):  
	''' Overlay additional shading onto a map

	Parameters
	----------
	dat : np.ndarrays with dimensions (y,x)
		Data to be plotted
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contourf except ``overlay``
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
	kwargs = __prepare_config(kwargs)

	dat = concat1(dat)

	def overlay(m, x, y, zorder, mask=None):
		# TODO: What about an additional colorbar for this data?
		__contourf_dat(m, x, y, dat, kwargs)

		return

	return overlay


def map_overlay_latlonbox(lon0, lon1, lat0, lat1, vertices=30, **kwargs):
	''' Overlay dots onto a map
	
	Parameters
	----------
	lon0, lat0 : float
		Location of the lower left corner
	lon1, lat1 : float
		Location of the upper right corner
	vertices : int
		How many interpolation points should be used to approximate the potentially
		curved sections of the latitude or longitude circles?

	Keyword arguments
	-----------------
	plot arguments : all contour
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''	
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


# that's it
