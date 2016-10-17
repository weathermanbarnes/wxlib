#!/usr/bin/env python3

import pickle
import os.path

from dynlib.shorthands import np, plt, fig
from dynlib.settings import get_active_context
from dynlib import proj

import dynlib.context.plot
conf = get_active_context()
u = ('u', 'u', 'U component of wind', 'm s**-1')
conf.register_variable([u, ], ['pv2000', ])

from create_ref_data import obfile, odfile, gfile


def gen_plot_setup_kwargs(name, **add_kwargs):
	def plot(dat, grid):
		plt.close('all')
		plotconf = conf.plot.merge(None, None, **add_kwargs)
		fig.setup(**plotconf)
		return fig.map(dat['pv']['u'][0,0,:,:], grid,
				_return=True, save='%s/%s.pdf' % (plotpath, name), show=False)
	
	return plot

def gen_plot_kwargs(name, **add_kwargs):
	def plot(dat, grid):
		plt.close('all')
		return fig.map(dat['pv']['u'][0,0,:,:], grid, 
				_return=True, save='%s/%s.pdf' % (plotpath, name), show=False,
				**add_kwargs)
	
	return plot

def gen_plot_contour_kwargs(name, **add_kwargs):
	def plot(dat, grid):
		plt.close('all')
		overlays = [
			fig.map_overlay_contour(dat['pv']['pt'][0,0,:,:], grid, **add_kwargs)
		]
		return fig.map(dat['pv']['u'][0,0,:,:], grid, m=proj.N_Atlantic, overlays=overlays,
				_return=True, save='%s/%s.pdf' % (plotpath, name), show=False)
	
	return plot


req_levs = set(['p', 'pv'])
# Currently untested: name, name_prefix, oroalpha, orocolor, oroscale, Zdata
# Currently no consequences: fig_dpi, hatches
plots = {
	# Contour(f) kwargs
	'default': gen_plot_kwargs('default'),
	'default_q': gen_plot_kwargs('default_q', q='u'),
	'default_plev': gen_plot_kwargs('default_plev', plev='pv2000'),
	'kwarg_alpha': gen_plot_kwargs('kwarg_alpha', alpha=0.5),
	'kwarg_cb_disabled': gen_plot_kwargs('kwarg_cb_disabled', cb_disabled=True),
	'kwarg_cb_expand_fig_fraction': gen_plot_kwargs('kwarg_cb_expand_fig_fraction', cb_expand_fig_fraction=0.5),
	'kwarg_cb_label': gen_plot_kwargs('kwarg_cb_label', cb_label='Wind speed [m/s]'),
	'kwarg_cb_orientation': gen_plot_kwargs('kwarg_cb_orientation', cb_orientation='vertical'),
	'kwarg_cb_tickspacing': gen_plot_kwargs('kwarg_cb_tickspacing', cb_tickspacing='uniform', scale=[0,30,50,60,80]),
	'kwarg_cb_shrink': gen_plot_kwargs('kwarg_cb_shrink', cb_shrink=0.5),
	'kwarg_cmap': gen_plot_kwargs('kwarg_cmap', cmap='viridis'),
	'kwarg_coastcolor': gen_plot_kwargs('kwarg_coastcolor', coastcolor='y'),
	'kwarg_colors': gen_plot_kwargs('kwarg_colors', scale=[0,30,50,70], colors=['r', 'y', 'g', 'c', 'b'], cmap=None),
	'kwarg_extend': gen_plot_kwargs('kwarg_extend', extend='max'),
	'kwarg_gridcolor': gen_plot_kwargs('kwarg_gridcolor', gridcolor='y'),
	'kwarg_hatches': gen_plot_kwargs('kwarg_hatches', scale=[0,30,50,70], colors=None, cmap=None, hatches=['.', '/', '\\', None, '*']),
	'kwarg_hook': gen_plot_kwargs('kwarg_hook', hook=lambda x: np.log(x)),
	'kwarg_m': gen_plot_kwargs('kwarg_m', m=proj.N_Atlantic),
	'kwarg_scale_colors': gen_plot_kwargs('kwarg_scale_colors', cmap=None, scale=[0,50], colors=['r', 'm', 'b']),
	'kwarg_scale_exceed_percentiles': gen_plot_kwargs('kwarg_scale_exceed_percentiles', scale_exceed_percentiles=(0.10,0.90)),
	'kwarg_scale_intervals': gen_plot_kwargs('kwarg_scale_intervals', scale_intervals=[3,]),
	'kwarg_scale_target_steps': gen_plot_kwargs('kwarg_scale_target_steps', scale_target_steps=20),
	'kwarg_scale_symmetric_zero': gen_plot_kwargs('kwarg_scale_symmetric_zero', scale_symmetric_zero=True),
	'kwarg_ticks': gen_plot_kwargs('kwarg_ticks', ticks=[0, 10, 50,]),
	'kwarg_ticklabels': gen_plot_kwargs('kwarg_ticklabels', ticks=[0, 10, 50,], ticklabels=['Zero', 'Little', 'Much']),
	'kwarg_tile': gen_plot_kwargs('kwarg_tile', tile=True),
	'kwarg_title': gen_plot_kwargs('kwarg_title', title='Custom title'),
	# Contour kwargs
	'oc_default': gen_plot_contour_kwargs('oc_default'),
	'oc_kwarg_alpha': gen_plot_contour_kwargs('oc_kwarg_alpha', alpha=0.5),
	'oc_kwarg_colors': gen_plot_contour_kwargs('oc_kwarg_colors', colors='b'),
	'oc_kwarg_contour_labels': gen_plot_contour_kwargs('oc_kwarg_contour_labels', contour_labels=True),
	'oc_kwarg_contour_labels_fontsize': gen_plot_contour_kwargs('oc_kwarg_contour_labels_fontsize', contour_labels=True, contour_labels_fontsize=20),
	'oc_kwarg_contour_labels_inline': gen_plot_contour_kwargs('oc_kwarg_contour_labels_inline', contour_labels=True, contour_labels_inline=False),
	'oc_kwarg_contour_labels_inline_spacing': gen_plot_contour_kwargs('oc_kwarg_contour_labels_inline_spacing', contour_labels=True, contour_labels_inline_spacing=10),
	'oc_kwarg_contour_labels_format': gen_plot_contour_kwargs('oc_kwarg_contour_labels_format', contour_labels=True, contour_labels_format='%.0f'),
	'oc_kwarg_linestyles': gen_plot_contour_kwargs('oc_kwarg_linestyles', linestyles='dashed'),
	'oc_kwarg_linewidths': gen_plot_contour_kwargs('oc_kwarg_linewidths', linewidths=0.1),
	# Setup kwargs
	'setup_default': gen_plot_setup_kwargs('setup_default'),
	'setup_kwarg_fig_dpi': gen_plot_setup_kwargs('setup_kwarg_fig_dpi', fig_dpi=300),
	'setup_kwarg_fig_size': gen_plot_setup_kwargs('setup_kwarg_fig_size', fig_size=(10,5)),
}
ofile = 'ref_plots.pickle'
plotpath = 'ref_plots'

if __name__ == '__main__':
	f = open(gfile, 'rb')
	grid = pickle.load(f)
	grid.oro = np.zeros(grid.x.shape)
	grid.cyclic_ew = True
	f.close()

	dat = {}
	
	# Get data to be plotted
	for pshort in req_levs:
		if not os.path.exists(obfile % pshort):
			raise IOError('Did not find cached base file for level `%s`.' % pshort)
		if not os.path.exists(odfile % pshort):
			raise IOError('Did not find cached derived file for level `%s`.' % pshort)
			
		print('Reading cached data for level %s' % pshort)

		dat[pshort] = {}

		f = np.load(obfile % pshort)
		dat[pshort].update(dict(f))
		f.close()
		f = np.load(odfile % pshort)
		dat[pshort].update(dict(f))
		f.close()

	# Do the plots
	imgs = {}
	size = 0.0
	cnt = 0
	for name, plot in plots.items():
		img = plot(dat, grid)
		size += img.getbuffer().nbytes
		cnt += 1
		imgs[name] = img
		
	# Save results as reference
	print('Saving plots (%d images, %.2fM)' % (cnt, size/1024**2) )
	f = open(ofile, 'wb')
	pickle.dump(imgs, f)
	f.close()

# C'est le fin
