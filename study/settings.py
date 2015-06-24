#!/usr/bin/env python
# -*- encoding: utf-8


''' Defines basic data structures for settings storage and convenient access

To be coded and written.
'''

import numpy as np
from copy import copy 

import lib.proj as proj
import lib.cm as cm

from collections import MutableMapping as mutmap


__context__ = set([])

def in_context(*contexts):
	for c in contexts:
		if c in __context__:
			return True
	
	return False

def def_context(context):
	__context__.add(context)




class default_dict(mutmap):
	def __init__(self, defaults):
		self._ = {}
		self._defaults = defaults

		mutmap.__init__(self)
		
	def __getitem__(self, key):
		if key in self._:
			return self._[key]
		elif key in self._defaults:
			return self._defaults[key]
		else:
			raise KeyError, key

	def __setitem__(self, key, value):
		if key in self._defaults:
			self._[key] = value
		else:
			raise KeyError, key

	def __delitem__(self, key):
		if key in self._:
			del self._[key]
		elif key in self._defaults:
			pass
		else:
			raise KeyError, key

	def __iter__(self):
		return self._defaults.__iter__()

	def __len__(self):
		return len(self._defaults)


class nd_default_dict(default_dict):
	""" n-dimensional dictionary with default values for valid keys in the last dimension

	The default value is independent from any dimension but the last, here called "first 
	dimensions".

	Valid combinations of keys for the first dimensions are specified using "tables". Each
	table is a list of sets containing valid keys for each of the first dimensions. Every 
	combination between these valid keys is assumed to be meaningful. As an example, an 
	array could contain a set of all available pressure levels and a set of all variables 
	available on pressure levels. If a variable is only available on one of the pressure
	levels, it would require an additional table, if one wants to avoid that configuration
	for that variable are available on all pressure levels.
	
	The set of allowed keys for each dimension is limited. Valid keys can be introduced 
	through the ``add_default()`` and ``add_table()`` functions.
	"""
	def __init__(self, first_table, defaults):
		self._fdims = [first_table, ] 	# List of tableSpecs, where each tableSpec is list of sets containing valid keys for each dimension
		self._ndims = len(first_table) + 1
		default_dict.__init__(self, defaults) 	# dictionary of valid keys in the last dimension and their default values
	
	def __check_key(self, key):
		if len(key) < self._ndims - 1:
			raise KeyError, 'Key must contain at least %d dimensions' % self._ndims-1
		elif len(key) > self._ndims:
			raise KeyError, 'Key can at most contain %d dimensions' % self._ndims
		
		found = False
		for table in self._fdims:
			for i, valid_skeys in zip(range(self._ndims-1), table):
				if not key[i] in valid_skeys:
					break
				elif i == self._ndims - 2:
					found = True
			if found == True:
				break
		if not found:
			raise KeyError, key
		
		return

	
	def __getitem__(self, key):
                self.__check_key(key)
		
		if len(key) == self._ndims:
			if key in self._: 
				return self._[key]
			elif key[-1] in self._defaults:
				return self._defaults[key[-1]]
			else:
				raise KeyError, key
		else: 
			ret = {}
			for skey in self._defaults:
				fkey = list(key) + [skey, ]
				ret[skey] = self[tuple(fkey)]

			return ret
	
	def __setitem__(self, key, value):
		self.__check_key(key)
		
		if len(key) == self._ndims:
			self._[key] = value
		else: 
			if not type(value) == dict:
				for lkey, lvalue in dict.items:
					fullkey = tuple(list(key)+[lkey,])
					if not lkey in self._defaults:
						raise KeyError, fullkey
					self._[fullkey] = lvalue
			else:
				raise ValueError, 'Dict required for multiple assignment.'

	def __delitem__(self, key):
		self.__check_key(key)
		
		if len(key) == self._ndims:
			if key in self._:
				del self._[key]
		else:
			for lkey in self._defaults:
				fullkey = tuple(list(key)+[lkey,])
				if fullkey in self._:
					del self._[fullkey]

	def add_table(self, table):
		if not len(table) == self._ndims - 1:
			raise ValueError, 'Table needs to have %d dimensions, got %d instead.' % (
					self._ndims - 1, len(table))
		self._fdims.append(table)

	def add_default(self, key, value):
		if key in self._defaults:
			raise KeyError, 'Default value for %s exists already' % str(key)

		self._defaults[key] = value



#class settings_dict(default_dict):
#	pass

class plot_settings_dict(nd_default_dict):
	def __init__(self, defaults):
		nd_default_dict.__init__(self, [], defaults)
		self._ndims = 3
		self._fdims = []


class settings_obj(default_dict):
	def _get(self, *keys):
		if len(keys) > 1:
			return self[keys[0]]._get(keys[1:])
		elif len(keys) == 1:
			return self[keys[0]]
		else:
			raise TypeError, '_get() requires at least one key'
	
	def _add_default(self, key, value):
		if key in self._defaults:
			raise KeyError, str(key) + ' already set'
		self._defaults[key] = value
	
	def __getattribute__(self, key):
		if key[0] == '_':
			return default_dict.__getattribute__(self, key)

		return self[key]

	def __setattr__(self, key, value):
		if key[0] == '_':
			return default_dict.__setattr__(self, key, value)
		
		self[key] = value
	
	def __delattr__(self, key):
		if key[0] == '_':
			return default_dict.__delattr__(self, key)
		
		del self[key]

	def register_variable(q, q_file=None, q_long='', q_unit=''):
		pass

	def register_level(plev):
		pass



PLOT_DEFAULTS = {
	'coastcolor': 'k', 
	'disable_cb': False, 
	'gridcolor': 'k',
	'hook': None,
	'm': proj.world, 
	'maskcolor': '0.25',
	'mark': None, 
	'oroalpha': 0.4, 
	'orocolor': 'k', 
	'oroscale': np.arange(1000,9000,1000),
	'overlays': [], 
	'plev': None, 
	'save': '', 
	'scale': 'auto', 
	'scale_exceed_percentiles': (0.01, 0.99),
	'scale_intervals': [1,2,3,5,10,],
	'scale_intervals_periodic': True,
	'scale_target_steps': 7,
	'scale_symmetric_zero': False,
	'show': True, 
	'ticks': None, 
	'ticklabels': [], 
	'title': '', 
	'Zdata': None,

	'lon': None, 
	'lat': None, 
	'alpha': 1.0, 
	'norm': None, 
	'vmin': None,
	'vmax': None, 
	'levels': None, 
	'origin': None, 
	'extent': None, 
}
PLOTF_DEFAULTS = copy(PLOT_DEFAULTS) 

PLOT_DEFAULTS.update({
	'cmap': None, 
	'colors': 'k', 
	'contour_labels': False,
	'contour_labels_fontsize': 12,
	'contour_labels_inline': True,
	'contour_labels_inline_spacing': 2,
	'contour_labels_format': '%1.1f',
	'extend': 'neither', 
	'linestyles': None, 
	'linewidths': 2.0, 
})
PLOTF_DEFAULTS.update({
	'cb_orientation': 'horizontal',
	'cmap': cm.gist_ncar, 
	'colors': None, 
	'extend': 'both', 

	'hatches': None 
})


conf = settings_obj({
	'q': [],
	'qi': [],
	'q_units': [],
	'q_long': [],
	'bins': [],
	'datapath': ['.', ],
	'opath': '.',
	'plotpath': '.',
	'file_std': '',
	'file_stat': '',
	'file_mstat': '',
	'years': [],
	'plevs': [],
	'ptlevs': [],
	'pvlevs': [],
	'zlevs': [],
	'mlevs': [],
	'plot': plot_settings_dict(PLOT_DEFAULTS),
	'plotf': plot_settings_dict(PLOTF_DEFAULTS),
})

# that's it
