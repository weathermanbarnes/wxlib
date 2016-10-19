#!/usr/bin/env python3

import pickle
import os.path

from lib.shorthands import np, plt, fig
from lib.settings import get_active_context
from lib import proj

import dynlib.context.plot
conf = get_active_context()
u = ('u', 'u', 'U component of wind', 'm s**-1')
conf.register_variable([u, ], ['pv2000', ])

from create_ref_data import path, obfile, odfile, gfile


ofile = 'ref_plots.pickle'
conf.plotpath = path+'/ref_plots'



def gen_plot_setup_kwargs(name, save, **add_kwargs):
	# Note: The save format here must be png, otherwise some of the comparisons will spuriously fail
	if save:
		save = '%s.png' % name
	
	def plot(dat, grid):
		plt.close('all')
		plotconf = conf.plot.merge(None, None, **add_kwargs)
		fig.setup(**plotconf)
		return fig.map(dat['pv']['u'][0,0,:,:], grid,
				_return=True, save=save, show=False)
	
	return plot

def gen_plot_kwargs(name, save, **add_kwargs):
	# Note: The save format here must be png, otherwise some of the comparisons will spuriously fail
	if save:
		save = '%s.png' % name

	def plot(dat, grid):
		plt.close('all')
		return fig.map(dat['pv']['u'][0,0,:,:], grid, 
				_return=True, save=save, show=False,
				**add_kwargs)
	
	return plot

def gen_plot_contour_kwargs(name, save, **add_kwargs):
	# Note: The save format here must be png, otherwise some of the comparisons will spuriously fail
	if save:
		save = '%s.png' % name

	def plot(dat, grid):
		plt.close('all')
		overlays = [
			fig.map_overlay_contour(dat['pv']['pt'][0,0,:,:], grid, **add_kwargs)
		]
		return fig.map(dat['pv']['u'][0,0,:,:], grid, m=proj.N_Atlantic, overlays=overlays,
				_return=True, save=save, show=False)
	
	return plot


req_levs = set(['p', 'pv'])
# Currently untested: name, name_prefix, oroalpha, orocolor, oroscale, Zdata
# Currently no consequences: fig_dpi, hatches
def make_plot_list(save):
	plots = {
		# Contour(f) kwargs
		'default': gen_plot_kwargs('default', save),
		'default_q': gen_plot_kwargs('default_q', save, q='u'),
		'default_plev': gen_plot_kwargs('default_plev', save, plev='pv2000'),
		'kwarg_alpha': gen_plot_kwargs('kwarg_alpha', save, alpha=0.5),
		'kwarg_cb_disabled': gen_plot_kwargs('kwarg_cb_disabled', save, cb_disabled=True),
		'kwarg_cb_expand_fig_fraction': gen_plot_kwargs('kwarg_cb_expand_fig_fraction', save, cb_expand_fig_fraction=0.5),
		'kwarg_cb_label': gen_plot_kwargs('kwarg_cb_label', save, cb_label='Wind speed [m/s]'),
		'kwarg_cb_orientation': gen_plot_kwargs('kwarg_cb_orientation', save, cb_orientation='vertical'),
		'kwarg_cb_tickspacing': gen_plot_kwargs('kwarg_cb_tickspacing', save, cb_tickspacing='uniform', scale=[0,30,50,60,80]),
		'kwarg_cb_shrink': gen_plot_kwargs('kwarg_cb_shrink', save, cb_shrink=0.5),
		'kwarg_cmap': gen_plot_kwargs('kwarg_cmap', save, cmap='viridis'),
		'kwarg_coastcolor': gen_plot_kwargs('kwarg_coastcolor', save, coastcolor='y'),
		'kwarg_colors': gen_plot_kwargs('kwarg_colors', save, scale=[0,30,50,70], colors=['r', 'y', 'g', 'c', 'b'], cmap=None),
		'kwarg_extend': gen_plot_kwargs('kwarg_extend', save, extend='max'),
		'kwarg_gridcolor': gen_plot_kwargs('kwarg_gridcolor', save, gridcolor='y'),
		'kwarg_hatches': gen_plot_kwargs('kwarg_hatches', save, scale=[0,30,50,70], colors=None, cmap=None, hatches=['.', '/', '\\', None, '*']),
		'kwarg_hook': gen_plot_kwargs('kwarg_hook', save, hook=lambda x: np.log(x)),
		'kwarg_m': gen_plot_kwargs('kwarg_m', save, m=proj.N_Atlantic),
		'kwarg_scale_colors': gen_plot_kwargs('kwarg_scale_colors', save, cmap=None, scale=[0,50], colors=['r', 'm', 'b']),
		'kwarg_scale_exceed_percentiles': gen_plot_kwargs('kwarg_scale_exceed_percentiles', save, scale_exceed_percentiles=(0.10,0.90)),
		'kwarg_scale_intervals': gen_plot_kwargs('kwarg_scale_intervals', save, scale_intervals=[3,]),
		'kwarg_scale_target_steps': gen_plot_kwargs('kwarg_scale_target_steps', save, scale_target_steps=20),
		'kwarg_scale_symmetric_zero': gen_plot_kwargs('kwarg_scale_symmetric_zero', save, scale_symmetric_zero=True),
		'kwarg_ticks': gen_plot_kwargs('kwarg_ticks', save, ticks=[0, 10, 50,]),
		'kwarg_ticklabels': gen_plot_kwargs('kwarg_ticklabels', save, ticks=[0, 10, 50,], ticklabels=['Zero', 'Little', 'Much']),
		'kwarg_tile': gen_plot_kwargs('kwarg_tile', save, tile=True),
		'kwarg_title': gen_plot_kwargs('kwarg_title', save, title='Custom title'),
		# Contour kwargs
		'oc_default': gen_plot_contour_kwargs('oc_default', save),
		'oc_kwarg_alpha': gen_plot_contour_kwargs('oc_kwarg_alpha', save, alpha=0.5),
		'oc_kwarg_colors': gen_plot_contour_kwargs('oc_kwarg_colors', save, colors='b'),
		'oc_kwarg_contour_labels': gen_plot_contour_kwargs('oc_kwarg_contour_labels', save, contour_labels=True),
		'oc_kwarg_contour_labels_fontsize': gen_plot_contour_kwargs('oc_kwarg_contour_labels_fontsize', save, contour_labels=True, contour_labels_fontsize=20),
		'oc_kwarg_contour_labels_inline': gen_plot_contour_kwargs('oc_kwarg_contour_labels_inline', save, contour_labels=True, contour_labels_inline=False),
		'oc_kwarg_contour_labels_inline_spacing': gen_plot_contour_kwargs('oc_kwarg_contour_labels_inline_spacing', save, contour_labels=True, contour_labels_inline_spacing=10),
		'oc_kwarg_contour_labels_format': gen_plot_contour_kwargs('oc_kwarg_contour_labels_format', save, contour_labels=True, contour_labels_format='%.0f'),
		'oc_kwarg_linestyles': gen_plot_contour_kwargs('oc_kwarg_linestyles', save, linestyles='dashed'),
		'oc_kwarg_linewidths': gen_plot_contour_kwargs('oc_kwarg_linewidths', save, linewidths=0.1),
		# Setup kwargs
		'setup_default': gen_plot_setup_kwargs('setup_default', save),
		'setup_kwarg_fig_dpi': gen_plot_setup_kwargs('setup_kwarg_fig_dpi', save, fig_dpi=300),
		'setup_kwarg_fig_size': gen_plot_setup_kwargs('setup_kwarg_fig_size', save, fig_size=(10,5)),
	}

	return plots

if __name__ == '__main__':
	f = open(path + '/' + gfile, 'rb')
	grid = pickle.load(f)
	grid.oro = np.zeros(grid.x.shape)
	grid.cyclic_ew = True
	f.close()

	dat = {}
	
	# Get data to be plotted
	for pshort in req_levs:
		if not os.path.exists(path + '/' + obfile % pshort):
			raise IOError('Did not find cached base file for level `%s`.' % pshort)
		if not os.path.exists(path + '/' + odfile % pshort):
			raise IOError('Did not find cached derived file for level `%s`.' % pshort)
			
		print('Reading cached data for level %s' % pshort)

		dat[pshort] = {}

		f = np.load(path + '/' + obfile % pshort)
		dat[pshort].update(dict(f))
		f.close()
		f = np.load(path + '/' + odfile % pshort)
		dat[pshort].update(dict(f))
		f.close()

	# Do the plots
	imgs = {}
	size = 0.0
	cnt = 0
	plots = make_plot_list(save=True)
	for name, plot in plots.items():
		img = plot(dat, grid)
		size += img.getbuffer().nbytes
		cnt += 1
		imgs[name] = img
		
	# Save results as reference
	print('Saving plots (%d images, %.2fM)' % (cnt, size/1024**2) )
	f = open(path + '/' + ofile, 'wb')
	pickle.dump(imgs, f)
	f.close()

# C'est le fin
