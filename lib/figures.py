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

from __future__ import absolute_import, unicode_literals, division, print_function
from six import string_types

import copy
import re
import types
from io import BytesIO
from PIL import Image

import numpy as np
import scipy.interpolate as intp

import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import LineCollection
import mpl_toolkits.basemap
basemap_version = mpl_toolkits.basemap.__version__

from . import settings as s

from .metio import metopen
from .utils import concat1, concat1lonlat, unflatten_lines, sect_gen_points
from .autoscale import autoscale


def _interpolate_for_section(static, dat, xlon, xlat):
	try:
		interp = intp.RectBivariateSpline(static.x[0,:], static.y[::-1,0], dat[::-1,:].T)
	# If interpolation fails, try again, but this time do not assume a structured mesh
	except TypeError:
		points = list(zip(static.x.flat, static.y.flat))
		interp = intp.LinearNDInterpolator(points, dat.flatten())
		points = list(zip(xlon, xlat))
		dati = interp(points)
	# If interpolation successful, get values for given points
	else:
		dati = interp.ev(xlon, xlat)
	
	return dati

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
	plev, q, kwargs = __prepare_config(kwargs)
	kwargs_map = copy.copy(kwargs)
	mask = __map_create_mask(static, kwargs)
	
	if 'map_axes' in kwargs:
		plt.sca(kwargs['map_axes'])
	else:
		plt.subplot(212)
	m, x, y, lon, lat = __map_setup(mask, static, kwargs)
	
	# Allow both homogenously increasing as well as homogeneously decreasing y axis coordinates
	if y[0,0] > y[1,0]:
		yslc = slice(None,None,-1)
	else:
		yslc = slice(None)
	
	# 1b. Create and interpolate points of the cross section
	xlon, xlat, xxy = sect_gen_points(sect, m, kwargs.pop('section_hor_resolution'))

	dati = np.empty((dat.shape[0], len(xlon),))
	for i in range(dat.shape[0]):
		dati[i] = _interpolate_for_section(static, dat[i,:,:], xlon, xlat)

	psi = _interpolate_for_section(static, ps, xlon, xlat)
	if not type(p) == type(None):
		pi = _interpolate_for_section(static, p, xlon, xlat)
	
	# 2a. Plot the inset map
	if not m or type(m) == type(plt):
		xx, xy = xlon, xlat
	else:
		xx, xy = m(xlon, xlat)
	if not type(datmap) == type(None): 
		datmap = __map_prepare_dat(datmap, mask, static, kwargs_map)
		__contourf_dat(m, x, y, datmap, q, kwargs_map)
	m.plot(xx, xy, 'k-', linewidth=4)
	
	logp = kwargs.get('logp', False)
	if logp:
		z = np.array([np.log(float(p_)) for p_ in static.z])
		kwargs['ticks'] = z
		kwargs['ticklabels'] = list(static.z)
	else:
		z = static.z
	
	# If pressure varies horizontally it needs to be interpolated as well
	if len(z.shape) > 1:
		zi = np.empty((z.shape[0], len(xlon),))
		for k in range(z.shape[0]):
			interp = intp.RectBivariateSpline(static.x[0,:], static.y[yslc,0], z[k,yslc,:].T)
			zi[k] = interp.ev(xlon, xlat)

		z = zi

	# 2b. Plot the actual cross section
	if 'section_axes' in kwargs:
		plt.sca(kwargs['section_axes'])
	else:
		plt.subplot(211)
	
	if kwargs.get('hook'):
		dati = kwargs.pop('hook')(dati)
	xxy = np.array(xxy)
	if len(xxy.shape) < len(z.shape):
		xxy_ = np.tile(xxy, (z.shape[0], 1))
	else:
		xxy_ = xxy
	__contourf_dat(plt, xxy_/1e3, z, dati, q, kwargs)
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
	plt.ylim(z.max(), z.min())
	plt.ylabel('Pressure [hPa]')
	plt.xlabel('Distance along section [km]')

	# 3. Finish off
	__decorate(m, xxy_, z, xlon, xlat, slice(None), plev, q, kwargs)
	return __output(plev, q, kwargs)


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

	if hasattr(plt, '__dynlib_latest_cs'):
		delattr(plt, '__dynlib_latest_cs')

	# 1. Prepare
	plev, q, kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)

	dat = __map_prepare_dat(dat, mask, static, kwargs)
	
	m, x, y, lon, lat = __map_setup(mask, static, kwargs)
	
	# 2. Plot the actual data
	__contourf_dat(m, x, y, dat, q, kwargs)

	# 3. Finish off
	__decorate(m, x, y, lon, lat, mask, plev, q, kwargs)
	return __output(plev, q, kwargs)


# Helper functions
def setup(**kwargs):
	''' Create a figure object, using some magic to automatically find a suitable size

	To find a suitable size, the following information is used

	 #. If fig_size is set to a tuple, the given fig_size will be used
	 #. If fig_size is set to a number, the aspect ratio of the map projection will be 
	    used to match the given size in square inches.
	 #. If fig_size is set to ``'auto'``, a fig_size of 10 in² is used as a default.
	
	For all of the above cases, the figure size is subsequently adapted: 

	 #. If a colorbar is to be added, the figure size is expanded accordingly to keep 
	    the size of the actual plot constant.
	 #. If a title is to be added, the figure size is expanded accordingly to keep 
	    the size of the actual plot constant.

	Keyword arguments
	-----------------
	plot arguments : all
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	matplotlib.figure.Figure
		Figure object with defined size and plot margins
	'''

	left = 0.01
	right = 0.99
	top = 0.99
	bottom = 0.01

	figsize = kwargs.pop('fig_size')
	dpi = kwargs.get('fig_dpi')

	# Default size: 32 in² at 150dpi gives about 20x20cm image at typical screen resolutions of 100 dpi
	if figsize == 'auto':
		figsize = 32.0 
	
	# Use the aspect ratio of the projection (if given, otherwise 1.5 is used) 
	# to find a good image size
	if type(figsize) == int or type(figsize) == float:
		aspect = kwargs.get('aspect', getattr(kwargs.get('m'), 'aspect', 1.5))
		height = round(np.sqrt(figsize/aspect),1)
		figsize = (aspect*height, height)

	elif not type(figsize) == tuple: 
		raise ValueError("fig_size must be either the string 'auto', a number or a tuple.")
	
	# Adapt figure size automatically if a color bar is added
	if not kwargs.get('cb_disable'): 
		expand = kwargs.get('cb_expand_fig_fraction')
		if kwargs.get('cb_orientation') == 'horizontal':
			figsize = (figsize[0], figsize[1]/(1.0-expand))
		else:
			figsize = (figsize[0]/(1.0-expand), figsize[1])
	
	# Adapt figure size and top margin automatically if a title is added
	if kwargs.get('title'):
		titlesize = 0.3 	# height of the title in inches

		height = figsize[1]
		top = height*top / (height + titlesize)
		figsize = (figsize[0], height + titlesize)

	fio = plt.figure(figsize=figsize, dpi=dpi)
	fio.subplotpars.update(left=left, right=right, top=top, bottom=bottom)

	return fio

def __prepare_config(kwargs):
	''' Make sure kwargs contains a complete contourf plot configuration, 
	filling undefined keys from the dynlib configuration.'''

	plev = kwargs.pop('plev', None)
	q = kwargs.pop('q', None)
	kwargs = s.conf.plotf.merge(plev, q, **kwargs)

	# cmap might be a function returing the cmap; if so generate it now!
	if 'cmap' in kwargs and type(kwargs['cmap']) == types.FunctionType:
		kwargs['cmap'] = kwargs['cmap']()
	
	return plev, q, kwargs

def __line_prepare_config(kwargs):
	''' Make sure kwargs contains a complete contour plot configuration, 
	filling undefined keys from the dynlib configuration.'''

	plev = kwargs.pop('plev', None)
	q = kwargs.pop('q', None)
	kwargs = s.conf.plot.merge(plev, q, **kwargs)
	
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
		if static.cyclic_ew:
			mask = concat1(mask)
		return mask

	plev = kwargs.pop('plev', None)
	datZ = kwargs.pop('Zdata', None)
	
	if plev and not type(datZ) == np.ndarray:
		f,datZ = metopen(s.conf.file_agg % {'agg': 'all', 'time': '%d-%d' % (s.conf.years[0],s.conf.years[-1]), 'plev': plev, 'q': 'Z'}, 'z', no_static=True)
		if f: f.close()
	if type(datZ) == np.ndarray:
		mask = datZ[:,:] < static.oro[:,:]
		if static.cyclic_ew:
			mask = concat1(mask)
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
		if isinstance(m, mpl_toolkits.basemap.Basemap):
			plt.gca().set_aspect('equal')
	if type(mask) == np.ndarray:
		m.contourf(x, y, mask, latlon=True, colors=kwargs.pop('maskcolor'))
		if isinstance(m, mpl_toolkits.basemap.Basemap):
			plt.gca().set_aspect('equal')
	
	return m, x, y, lon, lat

def __contourf_dat(m, x, y, dat, q, kwargs):
	''' Plot the actual data '''

	hatch = kwargs.pop('hatches')
	scale = kwargs.pop('scale')
	if isinstance(scale, string_types) and scale == 'auto':
		scale = autoscale(dat, **kwargs)
	
	if kwargs.get('tile'):
		if not 'edgecolors' in kwargs:
			kwargs['edgecolors'] = 'none'
		if not type(kwargs.get('colors')) == type(None) and type(kwargs.get('cmap')) == type(None):
			if not type(scale) == int:
				colors = list(kwargs.get('colors'))
				repeat = len(scale) // len(colors) + 1
				kwargs['cmap'] = matplotlib.colors.ListedColormap(colors*repeat)
			else:
				kwargs['cmap'] = matplotlib.colors.ListedColormap(kwargs.get('colors'))
		if not type(scale) == int:
			cmap = plt.get_cmap(kwargs['cmap'])
			#cmap.set_bad(alpha=0)
			#cmap.set_over(alpha=0)
			#cmap.set_under(alpha=0)
			kwargs['norm'] = matplotlib.colors.BoundaryNorm(scale, cmap.N)
			kwargs['cmap'] = cmap

		pkwargs = ['cmap', 'norm', 'vmin', 'vmax', 'edgecolors', 'alpha', 'clim', ]
		pkwargs = { key: kwargs.get(key, None) for key in pkwargs }
		datm = np.ma.masked_where(np.isnan(dat), dat)
		cs = m.pcolormesh(x, y, datm, latlon=True, zorder=1, **pkwargs)
	elif kwargs.get('tri'):
		cs = m.contourf(x.flatten(), y.flatten(), dat.flatten(), scale, latlon=True, zorder=1, **kwargs)
	else:
		cs = m.contourf(x, y, dat, scale, latlon=True, zorder=1, **kwargs)
	
	if isinstance(m, mpl_toolkits.basemap.Basemap):
		plt.gca().set_aspect('equal')

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
	
	if not kwargs.get('cb_disable'):
		plt.__dynlib_latest_cs = cs
		plt.__dynlib_latest_cs_kwargs = kwargs
		plt.__dynlib_latest_cs_q = q
	
	return

def __decorate(m, x, y, lon, lat, mask, plev, q, kwargs):
	''' Add "decorations": colorbar, legends, overlays and a title'''

	for overlay in kwargs.pop('overlays'):
		overlay(m, x,y, lon,lat, zorder=3, mask=mask)
	
	if hasattr(plt, '__dynlib_latest_cs'):
		cbcs = plt.__dynlib_latest_cs
		cbkwargs = plt.__dynlib_latest_cs_kwargs
		#cbq = plt.__dynlib_latest_cs_q <- unused?
		
		axes = cbkwargs.pop('cb_axes', None)
		orient = cbkwargs.pop('cb_orientation')
		spacing = cbkwargs.pop('cb_tickspacing')
		shrink = cbkwargs.pop('cb_shrink')
		expand = cbkwargs.pop('cb_expand_fig_fraction')
		padding = cbkwargs.pop('cb_expand_pad_fraction', 0.2)
		
		pad = expand * padding
		frac = expand * (1-padding)
	
		if axes:
			cb = plt.colorbar(cbcs, cax=axes, ticks=cbkwargs.pop('ticks'), orientation=orient)
		else:
			cb = plt.colorbar(cbcs, ticks=cbkwargs.pop('ticks'), orientation=orient, 
				shrink=shrink, pad=pad, fraction=frac, spacing=spacing)
		if cbkwargs.get('ticklabels'): 
			if not orient == 'vertical':
				cb.ax.set_xticklabels(cbkwargs.pop('ticklabels'))
			else:
				cb.ax.set_yticklabels(cbkwargs.pop('ticklabels'))
		if cbkwargs.get('cb_label'):
			if not orient == 'vertical':
				cb.ax.set_xlabel(cbkwargs.pop('cb_label'))
			else:
				cb.ax.set_ylabel(cbkwargs.pop('cb_label'))
	
	if kwargs.get('mark'):
		lons, lats = kwargs.pop('mark')
		mark_kwargs = kwargs.pop('mark_conf', dict(marker='o', facecolors=(1,1,1,0), 
				edgecolors='k', linewidths=3, size=484))
		size = mark_kwargs.pop('size')
		m.scatter(lons, lats, size, latlon=True, zorder=3, **mark_kwargs)
	
	if kwargs.get('title'):
		title = kwargs.pop('title')
		if title == 'auto':
			title = u'%s @ %s' % (s.conf.q_long.get(q, q), plev)
			if kwargs.get('name'):
				title += u' for %s' % kwargs.get('name')

		plt.title(title)
	
	return

def __output(plev, q, kwargs):
	''' Save and/or show the plot '''

	dpi = kwargs.pop('fig_dpi')
	if kwargs.get('save'):
		if hasattr(plt, '__dynlib_latest_cs_q'): 
			q = plt.__dynlib_latest_cs_q
		filename = kwargs.pop('save')
		if filename == 'auto':
			filename = '%s_%s_%s.%s' % (q, plev, __safename(kwargs.get('name', 'unnamed')), s.conf.plotformat)
		if kwargs.get('name_prefix'):
			filename = '%s_%s' % (kwargs.pop('name_prefix'), filename)
		
		# If png: Use adaptive palette to save space
		if filename[-3:] == 'png':
			imgstr = BytesIO()
			plt.savefig(imgstr, format='png', dpi=dpi)
			imgstr.seek(0)
			img = Image.open(imgstr)
			img_adaptive = img.convert('RGB').convert('P', palette=Image.ADAPTIVE)
			img_adaptive.save('%s/%s' % (s.conf.plotpath, filename), format='PNG')

		else:
			plt.savefig('%s/%s' % (s.conf.plotpath, filename), dpi=dpi)
	
	if kwargs.pop('show'):
		plt.show()
	
	if kwargs.pop('_return', False):
		imgstr = BytesIO()
		plt.savefig(imgstr, format='png', dpi=dpi)
		return imgstr
	
	return


def __safename(name):
	''' Make a plot name suitable and safe to be used as a file name segment '''

	name = name.lower()
	
	# Month names to numbers for sorting
	to_replace = {'jan': '01', 'feb': '02', 'mar': '03', 'apr': '04', 'mai': '05', 'jun': '06',
			'jul': '07', 'aug': '08', 'sep': '09', 'oct': '10', 'nov': '11', 'dec': '12'}
	# Replace some critical characters
	to_replace.update({u'æ': 'ae', u'ø': 'oe', u'å': 'aa', u'ä': 'ae', u'ö': 'oe', u'ü': 'ue', u'ß': 'ss'})

	for from_, to in to_replace.items():
		name = name.replace(from_, to)
	
	# Remove anything else that does not fit in the set {A-Z, a-z, 0-9, _, -, @}
	name = re.sub('([^A-Za-z0-9_\-@]+)', '_', name)

	return name
	


###############################################
# Overlays

def section_overlay_contour(dat, static, **kwargs):
	''' Overlay contours onto a section
	
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

	def overlay(m, xxy, z, xlon, xlat, zorder, mask=None):
		# Allow both homogenously increasing as well as homogeneously decreasing y axis coordinates
		if static.y[0,0] > static.y[1,0]:
			yslc = slice(None,None,-1)
		else:
			yslc = slice(None)

		dati = np.empty((dat.shape[0], len(xlon),))
		for i in range(dat.shape[0]):
			dati[i] = _interpolate_for_section(static, dat[i,:,:], xlon, xlat)
	
		if type(mask) == np.ndarray:
			dati[mask] = np.nan
		scale = kwargs.pop('scale')
		cs =  plt.contour(xxy/1e3, z, dati, scale, **kwargs)

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


def map_overlay_contour(dat, static, **kwargs):  
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

	if static.cyclic_ew:
		x_, y_ = concat1lonlat(static.x, static.y)
	else:
		x_, y_ = static.x, static.y
	
	kwargs = __line_prepare_config(kwargs)
	owngrid = kwargs.pop('owngrid', False)
	
	def overlay(m, x, y, lon, lat, zorder, mask=None):
		if owngrid:
			mask = __map_create_mask(static, kwargs)
			dat_ = __map_prepare_dat(dat, mask, static, kwargs)
			m, x, y, lon, lat = __map_setup(mask, static, kwargs)
		else:
			dat_ = __map_prepare_dat(dat, mask, static, kwargs)

		if type(mask) == np.ndarray:
			dat_[mask] = np.nan
		scale = kwargs.pop('scale')
		if isinstance(scale, string_types) and scale == 'auto':
			scale = autoscale(dat_, **kwargs)
		cs =  m.contour(x_, y_, dat_, scale, latlon=True, **kwargs)
		if isinstance(m, mpl_toolkits.basemap.Basemap):
			plt.gca().set_aspect('equal')

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

	lns = unflatten_lines(lines, loff, static)

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		# TODO: Convert to latlon=True system
		for ln in lns:
			xfr, yfr = m(ln[:,0], ln[:,1])
			m.plot(xfr, yfr, kwargs['linecolor'], linewidth=kwargs.get('linewidth', 2), 
					alpha=kwargs.get('alpha', 1))

		return

	return overlay


def map_overlay_clines(lines, cdat, loff, static, **kwargs):
	''' Overlay colorful lines onto a map

	In comparison with ``map_overlay_lines``, this function takes an additional
	cdat argument providing the data basis for coloring the lines using a color map.
	
	Parameters
	----------
	lines : np.ndarray with dimensions (pointindex,infotype)
		Line position array
	cdat : np.ndarray with dimensions (pointindex)
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
	
	lns = unflatten_lines(lines, loff, static, 
		convert_grididx=kwargs.get('convert_grididx2lonlat', True))
	
	scale = kwargs.get('scale', None)
	if type(scale) in [list, tuple, np.ndarray,]:
		norm = plt.Normalize(scale[0], scale[-1])
	else:
		norm = plt.Normalize()
		norm.autoscale(cdat)

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		for lidx, ln in zip(range(len(lns)), lns):
			xfr, yfr = m(ln[:,0], ln[:,1])

			points = np.array([xfr, yfr]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)

			lc = LineCollection(segments, array=cdat[loff[lidx]:loff[lidx+1]],
					cmap=kwargs.get('cmap'), norm=norm,
					linewidth=kwargs.get('linewidth'), alpha=kwargs.get('alpha'),
			)
			plt.gca().add_collection(lc)

			if not kwargs.get('cb_disable'):
				plt.__dynlib_latest_cs = lc
				plt.__dynlib_latest_cs_kwargs = kwargs
				plt.__dynlib_latest_cs_q = kwargs.get('q')
		
		return

	return overlay



def map_overlay_dots(lons, lats, static, **kwargs):  
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

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		# TODO: Convert to latlon=True system
		xfr, yfr = m(lons, lats)
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

	plev, q, kwargs = __prepare_config(kwargs)

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		dat_ = __map_prepare_dat(dat, mask, static, kwargs)

		__contourf_dat(m, x, y, dat_, q, kwargs)

		return

	return overlay

def map_overlay_barbs(u, v, static, **kwargs):  
	''' Overlay wind barbs onto a map

	For rotated grids (as far as supported by gridlib), vectors will be 
	rotated back automatically.

	Parameters
	----------
	u : np.ndarrays with dimensions (y,x)
		x-component of the vector to be plotted.
	v : np.ndarrays with dimensions (y,x)
		y-component of the vector to be plotted.
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contour except ``overlay``
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
	
	kwargs = __line_prepare_config(kwargs)

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		def rotate_vector(u, v, lon, lat, kwargs):
			interval = kwargs.pop('vector_space_interval', 15)
			slc = (slice(interval//2,None,interval), slice(interval//2,None,interval))
			return m.rotate_vector(u_[slc], v_[slc], lon[slc], lat[slc], returnxy=True)


		u_ = __map_prepare_dat(u, mask, static, kwargs)
		v_ = __map_prepare_dat(v, mask, static, kwargs)

		# Respect rotated coordinate systems (otherweise returned unchanged)
		u_, v_ = static.unrotate_vector(u_, v_)
		
		# TODO: Why does transform_vector lead to errors for iveret!?
		if not kwargs.get('vector_disable_interpolation', False):
			try:
				if lat[0,0] > lat[-1,0]:
					ut,vt, xt,yt = m.transform_vector(u_[::-1,:],v_[::-1,:],lon[0,:],lat[::-1,0], 30, 20, returnxy=True)
				else:
					ut,vt, xt,yt = m.transform_vector(u_,v_,lon[0,:],lat[:,0], 30, 20, returnxy=True)

			except (AttributeError, ValueError):
				ut,vt, xt,yt = rotate_vector(u_, v_, lon, lat, kwargs)
		else:
			ut,vt, xt,yt = rotate_vector(u_, v_, lon, lat, kwargs)
				
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5, zorder=3)
	
	return overlay

def map_overlay_quiver(u, v, static, **kwargs):
	''' Overlay quiver vectors onto a map

	For rotated grids (as far as supported by gridlib), vectors will be 
	rotated back automatically.

	Parameters
	----------
	u : np.ndarrays with dimensions (y,x)
		x-component of the vector to be plotted.
	v : np.ndarrays with dimensions (y,x)
		y-component of the vector to be plotted.
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contour except ``overlay``
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''
	
	kwargs = __line_prepare_config(kwargs)

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		u_ = __map_prepare_dat(u, mask, static, kwargs)
		v_ = __map_prepare_dat(v, mask, static, kwargs)
		
		# Respect rotated coordinate systems (otherweise returned unchanged)
		u_, v_ = static.unrotate_vector(u_, v_)

		try:
			ut,vt, xt,yt = m.transform_vector(u_[::-1,:],v_[::-1,:],lon[0,:],lat[::-1,0], 30, 20, returnxy=True)
		except (AttributeError, ValueError):
			interval = kwargs.pop('vector_space_interval', 15)
			slc = (slice(interval//2,None,interval), slice(interval//2,None,interval))
			ut,vt, xt,yt = m.rotate_vector(u_[slc], v_[slc], lon[slc], lat[slc], returnxy=True)
		
		m.quiver(xt, yt, ut, vt, zorder=3, scale=kwargs.pop('quiver_length', None), scale_units='width',
				linewidths=kwargs.pop('linewidths', None))
	
	return overlay

def map_overlay_dilatation(defabs, defang, static, **kwargs):
	''' Overlay dilatation axes onto a map

	For rotated grids (as far as supported by gridlib), orientations will be 
	rotated back automatically.

	Parameters
	----------
	defabs : np.ndarrays with dimensions (y,x)
		Magnitude of the dilatation
	defang : np.ndarrays with dimensions (y,x)
		Orientation of the dilatation
	static : gridlib.grid
		Meta information about the data array, like the grid definition
	
	Keyword arguments
	-----------------
	plot arguments : all contour except ``overlay``
		For a list of valid arguments refer to :ref:`plot configuration`.
	
	Returns
	-------
	function
		Overlay as a callable function
	'''

	kwargs = __line_prepare_config(kwargs)

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		defdex = np.cos(defang[:,:]) *defabs
		defdey = np.sin(defang[:,:]) *defabs
	
		defdex = __map_prepare_dat(defdex, mask, static, copy.copy(kwargs))
		defdey = __map_prepare_dat(defdey, mask, static, kwargs)

		# Respect rotated coordinate systems (otherweise returned unchanged)
		defdex, defdey = static.unrotate_vector(defdex, defdey)

		# TODO: Better autoscaling, make scaling configurable!

		try:
			Nvecx = 27
			Nvecy = 18
			ut,vt,xt,yt = m.transform_vector(defdex[::-1,:],defdey[::-1,:],lon[0,:],lat[::-1,0], Nvecx, Nvecy, returnxy=True)
			qscale = 420
		except AttributeError:
			interval = kwargs.pop('vector_space_interval', 15)
			slc = (slice(interval//2,None,interval), slice(interval//2,None,interval))
			ut,vt,xt,yt = defdex[slc],defdey[slc],x[slc],y[slc]
			qscale = 72
		except ValueError:
			interval = kwargs.pop('vector_space_interval', 15)
			slc = (slice(interval//2,None,interval), slice(interval//2,None,interval))
			ut,vt,xt,yt = m.rotate_vector(defdex[slc],defdey[slc],lon[slc],lat[slc], returnxy=True)
			qscale = 840

		m.quiver(xt, yt, ut, vt, zorder=4, scale=qscale, alpha=0.85, 
				linewidths=kwargs.get('linewidths', None))
		m.quiver(xt, yt, -ut, -vt, zorder=4, scale=qscale, alpha=0.85, 
				linewidths=kwargs.pop('linewidths', None))

	return overlay

def map_overlay_latlonbox(lon0, lon1, lat0, lat1, vertices=30, **kwargs):
	''' Overlay dots onto a map
	
	Parameters
	----------
	lon0, lon1 : float or None
		Longitude of the west and east boundary
	lat0, lat1 : float or None
		Latitude of the north and south boundary
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

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		if lat0:
			lat0_ = lat0
		else:
			lat0_ = -90.0
		if lat1:
			lat1_ = lat1
		else:
			lat1_ = 90.0

		# Western boundary
		if lon0:
			m.plot(np.ones((vertices,))*lon0, np.linspace(lat0_,lat1_,vertices), latlon=True, **kwargs)
			lon0_ = lon0
		else: 
			lon0_ = -180.0
		# Eastern boundary
		if lon1:
			m.plot(np.ones((vertices,))*lon1, np.linspace(lat0_,lat1_,vertices), latlon=True, **kwargs)
			lon1_ = lon1
		else:
			lon1_ = 180.0

		# Southern boundary
		if lat0:
			m.plot(np.linspace(lon0_,lon1_,vertices), np.ones((vertices,))*lat0, latlon=True, **kwargs)
		# Northern boundary
		if lat1:
			m.plot(np.linspace(lon0_,lon1_,vertices), np.ones((vertices,))*lat1, latlon=True, **kwargs)

		return

	return overlay


# that's it
