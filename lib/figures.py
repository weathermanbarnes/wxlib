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
import re
import types
import cStringIO
from PIL import Image

import numpy as np
import scipy.interpolate as intp

import matplotlib.pyplot as plt
import matplotlib.colors
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
	plev, q, kwargs = __prepare_config(kwargs)
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
	__decorate(m, xxy/1e3, z, None, None, slice(None), plev, q, kwargs)
	__output(plev, q, kwargs)

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
	plev, q, kwargs = __prepare_config(kwargs)
	mask = __map_create_mask(static, kwargs)

	dat = __map_prepare_dat(dat, mask, static, kwargs)
	
	m, x, y, lon, lat = __map_setup(mask, static, kwargs)
	
	# 2. Plot the actual data
	__contourf_dat(m, x, y, dat, kwargs)

	# 3. Finish off
	__decorate(m, x, y, lon, lat, mask, plev, q, kwargs)
	__output(plev, q, kwargs)

	return


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
		aspect = getattr(kwargs.get('m'),'aspect', 1.5)
		height = round(np.sqrt(figsize/aspect),1)
		figsize = (aspect*height, height)

	elif not type(figsize) == tuple: 
		raise ValueError, "fig_size must be either the string 'auto', a number or a tuple."

	# Adapt figure size automatically if a color bar is added
	if not kwargs.get('cb_disabled'): 
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
	kwargs = conf.plotf.merge(plev, q, **kwargs)

	# cmap might be a function returing the cmap; if so generate it now!
	if 'cmap' in kwargs and type(kwargs['cmap']) == types.FunctionType:
		kwargs['cmap'] = kwargs['cmap']()
	
	return plev, q, kwargs

def __line_prepare_config(kwargs):
	''' Make sure kwargs contains a complete contour plot configuration, 
	filling undefined keys from the dynlib configuration.'''

	plev = kwargs.pop('plev', None)
	q = kwargs.pop('q', None)
	kwargs = conf.plot.merge(plev, q, **kwargs)
	
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
		f,datZ = metopen(conf.file_agg % {'agg': 'all', 'time': '%d-%d' % (conf.years[0],conf.years[-1]), 'plev': plev, 'q': 'Z'}, 'z', no_static=True)
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

def __decorate(m, x, y, lon, lat, mask, plev, q, kwargs):
	''' Add "decorations": colorbar, legends, overlays and a title'''

	if not kwargs.pop('cb_disable'):
		orient = kwargs.pop('cb_orientation')
		spacing = kwargs.pop('cb_tickspacing')
		shrink = kwargs.pop('cb_shrink')
		expand = kwargs.pop('cb_expand_fig_fraction')
		
		pad = expand/5
		frac = 4*expand/5

		cb = plt.colorbar(ticks=kwargs.pop('ticks'), orientation=orient, 
				shrink=shrink, pad=pad, fraction=frac, spacing=spacing)
		if kwargs.get('ticklabels'): 
			if not orient == 'vertical':
				cb.ax.set_xticklabels(kwargs.pop('ticklabels'))
			else:
				cb.ax.set_yticklabels(kwargs.pop('ticklabels'))
		if kwargs.get('cb_label'):
			if not orient == 'vertical':
				cb.ax.set_xlabel(kwargs.pop('cb_label'))
			else:
				cb.ax.set_ylabel(kwargs.pop('cb_label'))
	
	#legend_labels = kwargs.pop('legend_labels', None)
	#if legend_labels:
	#	plt.legend([], legend_labels)
	
	for overlay in kwargs.pop('overlays'):
		overlay(m, x,y, lon,lat, zorder=3, mask=mask)
	
	if kwargs.get('mark'):
		yidx, xidx = kwargs.pop('mark')
		m.scatter(x[yidx,xidx], y[yidx,xidx], 484, latlon=True, marker='o', facecolors=(1,1,1,0), 
				edgecolors='k', linewidths=3, zorder=3)
	
	if kwargs.get('title'):
		title = kwargs.pop('title')
		if title == 'auto':
			title = u'%s @ %s' % (conf.q_long.get(q, q), plev)
			if kwargs.get('name'):
				title += u' for %s' % kwargs.get('name')

		plt.title(title)
	
	return

def __output(plev, q, kwargs):
	''' Save and/or show the plot '''

	if kwargs.get('save'):
		filename = kwargs.pop('save')
		if filename == 'auto':
			filename = '%s_%s_%s.%s' % (q, plev, __safename(kwargs.get('name', 'unnamed')), conf.plotformat)
		if kwargs.get('name_prefix'):
			filename = '%s_%s' % (kwargs.pop('name_prefix'), filename)
		
		# If png: Use adaptive palette to save space
		if filename[-3:] == 'png':
			imgstr = cStringIO.StringIO()
			plt.savefig(imgstr, format='png', dpi=kwargs.pop('fig_dpi'))
			imgstr.seek(0)
			img = Image.open(imgstr)
			img_adaptive = img.convert('RGB').convert('P', palette=Image.ADAPTIVE)
			img_adaptive.save('%s/%s' % (conf.plotpath, filename), format='PNG')

		else:
			plt.savefig('%s/%s' % (conf.plotpath, filename), dpi=kwargs.pop('fig_dpi'))
	
	if kwargs.pop('show'):
		plt.show()
	
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

# TODO: Should this one take the static object as an argument? If only for consistency in the API.
def section_overlay_contour(dat, sect, **kwargs):
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

	kwargs = __line_prepare_config(kwargs)
	
	def overlay(m, x, y, lon, lat, zorder, mask=None):
		dat_ = __map_prepare_dat(dat, mask, static, kwargs)

		if type(mask) == np.ndarray:
			dat_[mask] = np.nan
		scale = kwargs.pop('scale')
		if scale == 'auto':
			scale = autoscale(dat_, **kwargs)
		cs =  m.contour(x, y, dat_, scale, latlon=True, **kwargs)

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

	def overlay(m, x, y, lon, lat, zorder, mask=None):
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

	def overlay(m, x, y, lon, lat, zorder, mask=None):
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

	def overlay(m, x, y, lon, lat, zorder, mask=None):
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

	plev, q, kwargs = __prepare_config(kwargs)

	def overlay(m, x, y, lon, lat, zorder, mask=None):
		dat_ = __map_prepare_dat(dat, mask, static, kwargs)

		# TODO: What about an additional colorbar for this data?
		__contourf_dat(m, x, y, dat_, kwargs)

		return

	return overlay

def map_overlay_barbs(u, v, static, **kwargs):  
	''' Overlay wind barbs onto a map

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

		try:
			ut,vt, xt,yt = m.transform_vector(u_[::-1,:],v_[::-1,:],lon[0,:],lat[::-1,0], 30, 20, returnxy=True)
		except ValueError:
			interval = kwargs.pop('vector_space_interval', 15)
			slc = (slice(interval/2,None,interval), slice(interval/2,None,interval))
			ut,vt, xt,yt = m.rotate_vector(u_[slc], v_[slc], lon[slc], lat[slc], returnxy=True)
		
		m.barbs(xt, yt, ut, vt, length=6, linewidth=0.5, zorder=3)
	
	return overlay

def map_overlay_quiver(u, v, static, **kwargs):  
	''' Overlay quiver vectors onto a map

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

		try:
			ut,vt, xt,yt = m.transform_vector(u_[::-1,:],v_[::-1,:],lon[0,:],lat[::-1,0], 30, 20, returnxy=True)
		except ValueError:
			interval = kwargs.pop('vector_space_interval', 15)
			slc = (slice(interval/2,None,interval), slice(interval/2,None,interval))
			ut,vt, xt,yt = m.rotate_vector(u_[slc], v_[slc], lon[slc], lat[slc], returnxy=True)
		
		m.quiver(xt, yt, ut, vt, zorder=3, scale=kwargs.pop('quiver_length', None), scale_units='width')
	
	return overlay

def map_overlay_dilatation(defabs, defang, static, **kwargs):
	''' Overlay dilatation axes onto a map

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

		Nvecx = 27
		Nvecy = 18
		
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

	def overlay(m, x, y, lon, lat, zorder, mask=None):
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
