#!/usr/bin/env python
# -*- encoding: utf-8


''' Defines basic data structures for settings storage and convenient access

Dynlib settings are in essence a (rather large) set of key-value pairs, which 
could be well represented by the python built-in ``dict`` object. However, the 
pure ``dict`` object has some drawbacks in this context:

 #. ``dict`` objects cannot represent default values for keys. Whenever a key in
    ``dict`` is overwritten, the original content is lost.
 #. ``dict`` cannot represent interdependencies between different keys.
 #. ``dict`` cannot represent several variants of itself. As an example, the 
    plot configuration for different variables will be largely identical, but only
    differ in a few configuration keys. When represented by pure ``dict`` objects, 
    the plot configuration for each variable would need to be full independent 
    from each other. This indepence complicates would makes it tedious and 
    error-prone to apply a customised configuration to all variables.

For these reasons, this module introduces:

 #. ``default_dict`` objects: A variant of ``dict`` with a pre-defined set of keys 
    and default values for each.
 #. ``nd_default_dict`` objects: Derived from ``default_dict``, it adds an abitrary
    number of dimensions, accessible by ``dat[key1,key2,...,keyN]``. Valid keys for
    each dimension are prescribed, and default values must be provided for each key
    in the last dimesion. Default values are identical for all ``dat[::,keyN]``, 
    and hence independent from anything but the last dimension.
 #. ``plot_settings_dict``: A variant of the ``nd_default_dict`` used for storing
    plot configuration. The number of dimensions is fixed to three. They are, in order: 
    vertical level, variable and plot configuration key. 
 #. ``settings_obj``: A variant of the ``default_dict`` where access is also possible
    through the attribute syntax dat.key. This class is the basis for the conf 
    object, the configuration root. Within, ``conf.plot`` and ``conf.plotf`` are 
    ``plot_settings_dict`` objects.
'''

import numpy as np
from copy import copy 

import proj
import cm

from collections import MutableMapping as mutmap


__context__ = set([])

def in_context(*contexts):
	''' Check if any of the given contexts is active
	
	Parameters
	----------
	context1 - contextN : str
	    Contexts to be checked
	
	Returns
	-------
	bool
	    ``True`` if any of the given contexts are active.
	'''

	for c in contexts:
		if c in __context__:
			return True
	
	return False

def def_context(context):
	''' Define a new context

	Parameters
	----------
	context : str
	    Name of the context to be defined
	'''

	__context__.add(context)


class default_dict(mutmap):
	''' Dictionary with a predefined set of valid keys and default values for each 

	Non-default values are stored in ``self._``, the defaults in ``self._defaults``.
	
	Parameters
	----------
	defaults : dict
	    Definition of valid keys and their default values
	
	Notes
	-----
	 #. ToDo: Re-Implement mutex groups
	'''

	def __init__(self, defaults):
		self._ = {}
		self._defaults = defaults

		mutmap.__init__(self)
		
	def __getitem__(self, key):
		''' Implements the access syntax ``dat[key]`` '''

		if key in self._:
			return self._[key]
		elif key in self._defaults:
			return self._defaults[key]
		else:
			raise KeyError, key

	def __setitem__(self, key, value):
		''' Implements the assignment syntax ``dat[key] = value`` '''

		if key in self._defaults:
			self._[key] = value
		else:
			raise KeyError, key

	def __delitem__(self, key):
		''' Implements the reset syntax ``del dat[key]`` '''

		if key in self._:
			del self._[key]
		elif key in self._defaults:
			pass
		else:
			raise KeyError, key

	def __iter__(self):
		''' Implement interation over the dictionary, e.g. the syntax ``for key in dat:`` '''

		return self._defaults.__iter__()

	def __len__(self):
		''' Implements the length query, allowing among others to use ``len(dat)`` '''

		return len(self._defaults)


class nd_default_dict(default_dict):
	''' n-dimensional dictionary with default values for valid keys in the last dimension

	The default value is independent from any dimension but the last. All dimensions 
	before the last are subsequently called "first dimensions".

	Valid combinations of keys for the first dimensions are specified using "tables". Each
	table is a list of sets containing valid keys for each of the first dimensions. Every 
	combination between these valid keys is assumed to be meaningful. As an example, an 
	array could contain a set of all available pressure levels and a set of all variables 
	available on pressure levels. If a variable is only available on one of the pressure
	levels, it would require an additional table, if one wants to avoid that configuration
	for that variable are available on all pressure levels.
	
	The set of allowed keys for each dimension is limited. Valid keys can be introduced 
	through the ``add_default()`` and ``add_table()`` functions.

	Parameters
	----------
	first_table : list of sets
	    Sets of valid keys for each of the first dimensions
	defaults : dict
	    Definition of valid keys and their default values
	
	Notes
	----
	 #. ToDo: Implement mutex groups
	'''

	def __init__(self, first_table, defaults):
		self._fdims = [first_table, ]
		self._ndims = len(first_table) + 1
		default_dict.__init__(self, defaults) 	# dictionary of valid keys in the last dimension and their default values
	
	def __check_key(self, key, allow_blank=False):
		''' Check if a given query key is valid '''

		if len(key) < self._ndims - 1:
			raise KeyError, 'Key must contain at least %d dimensions' % self._ndims-1
		elif len(key) > self._ndims:
			raise KeyError, 'Key can at most contain %d dimensions' % self._ndims
		
		found = False
		for table in self._fdims:
			for i, valid_skeys in zip(range(self._ndims-1), table):
				if key[i] == slice(None):
					if allow_blank:
						if i == self._ndims - 2:
							found = True
					else:
						raise KeyError, 'Only multiple assignment supported.'
				elif not key[i] in valid_skeys:
					break
				elif i == self._ndims - 2:
					found = True
			if found == True:
				break
		if not found:
			raise KeyError, key
		
		return

	def __resolve_keys(self, key):
		''' Find all keys that match a given query key for the first dimensions

		If the query key does not contain blanks, ``__resolve_keys`` will return
		a list with only that one key.
		'''

		# If key itself is not valid, return empty list
		try:
			self.__check_key(key, allow_blank=True)
		except KeyError:
			return []
		
		# Solve blank dimensions in the key
		keys = []
		blank_dim = False
		for n, dim in zip(range(self._ndims-1), key):
			if dim == slice(None):
				blank_dim = True
				
				# Find key fragments to substitute the blank
				valid_skeys = []
				for fdim in self._fdims:
					nofit = False
					for m in range(n):
						if key[m] == slice(None):
							continue
						elif not key[m] in fdim[m]:
							nofit = True
							break
					if not nofit:
						valid_skeys.extend(fdim[n])
				
				# Recursively solve for later dimensions
				for skey in valid_skeys:
					_key = key[:n] + (skey,) + key[n+1:]
					keys.extend(self.__resolve_keys(_key))
				
				# If there are more blank dimensions, 
				# they will be resolved in the recursive calls, not here!
				break
		
		# Alternatively, if there is no blank dimension, return the given key itself
		if not blank_dim:
			keys.append(key)

		return keys

	
	def __getitem__(self, key):
		''' Implement access via the ``dat[key1,key2,...,keyN]`` syntax '''

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
				fkey = key + (skey, )
				ret[skey] = self[fkey]

			return ret
	
	def __setitem__(self, key, value):
		''' Implement assignment via the ``dat[key1,key2,...,keyN] = value`` syntax '''

		self.__check_key(key, allow_blank=True)
		keys = self.__resolve_keys(key[:self._ndims-1])

		if len(key) == self._ndims:
			keys = [_key+(key[-1],) for _key in keys]
		
		for key in keys:
			if len(key) == self._ndims:
				self._[key] = value
			else: 
				if not type(value) == dict:
					for lkey, lvalue in dict.items:
						fullkey = key + (lkey,)
						if not lkey in self._defaults:
							raise KeyError, fullkey
						self._[fullkey] = lvalue
				else:
					raise TypeError, 'Dict required for multiple assignment.'

	def __delitem__(self, key):
		''' Implement reset via the ``del dat[key1,key2,...,keyN]`` syntax '''

		self.__check_key(key)
		
		if len(key) == self._ndims:
			if key in self._:
				del self._[key]
		else:
			for lkey in self._defaults:
				fullkey = key + (lkey,)
				if fullkey in self._:
					del self._[fullkey]

	def add_table(self, table):
		''' Add new combinations for the first dimensions '''

		if not len(table) == self._ndims - 1:
			raise ValueError, 'Table needs to have %d dimensions, got %d instead.' % (
					self._ndims - 1, len(table))
		self._fdims.append(table)

	def add_default(self, key, value):
		''' Add a new key and its default value for the last dimension '''

		if key in self._defaults:
			raise KeyError, 'Default value for %s exists already' % str(key)

		self._defaults[key] = value



#class settings_dict(default_dict):
#	pass

class plot_settings_dict(nd_default_dict):
	''' A version of the nd_default_dict, where the number of dimensions is fixed to 3

	The dimensions are, in order: vertical level, variable and plot configuration key.

	In contrast to the nd_default_dict, plot_settings_dict does not take an initial table as
	an argument, as the dimensions are prescribed and no (plev,q)-table is more equal than
	the others.

	Parameters
	----------
	defaults : dict
	    Definition of valid keys and their default values
	'''
	def __init__(self, defaults):
		nd_default_dict.__init__(self, [], defaults)
		self._ndims = 3
		self._fdims = []


class settings_obj(default_dict):
	''' Another interface to a dict using the attribute syntax for access 

	The object provides the primary access to dynlib configuration via its 
	attributes, e.g. ``conf.datapath``. As the configuration keys must be 
	valid attribute names, keys are restricted to strings containing valid
	python variable names. They can for example not begin with a number.

	Even though access through attributes is the prefered way of accessing
	configuration, the dictionary syntax ``conf['datapath']`` is not 
	disabled and provided identical functionality.

	It also provides a mechanism to define new variables/vertical levels 
	through the ``register_variable`` method.

	Parameters
	----------
	defaults : dict
	    Definition of valid keys and their default values
	'''
	
	# Configuration keys that have a special meaning for the functionality of this class:
	# They are required for registering variables
	_Q = 'q'
	_Q_FILE = 'qf'
	_Q_LONG = 'q_long'
	_Q_UNITS = 'q_units'
	_Q_BINS = 'q_bins'
	
	_PLOT = 'plot'
	_PLOTF = 'plotf'

	def _get(self, *keys):
		''' A third way to get a configuration item, using the syntax ``conf._get('datapath')`` '''

		if len(keys) > 1:
			return self[keys[0]]._get(keys[1:])
		elif len(keys) == 1:
			return self[keys[0]]
		else:
			raise TypeError, '_get() requires at least one key'
	
	def _add_default(self, key, value):
		''' Add a new configuration item and its default value '''

		if key in self._defaults:
			raise KeyError, str(key) + ' already set'
		self._defaults[key] = value
	
	def __getattribute__(self, key):
		''' Implement the attribute access syntax ``conf.datapath`` '''

		if key[0] == '_' or key == 'register_variable':
			return default_dict.__getattribute__(self, key)

		return self[key]

	def __setattr__(self, key, value):
		''' Implement the attribute assignment syntax ``conf.datapath = value`` '''

		if key[0] == '_' or key == 'register_variable':
			return default_dict.__setattr__(self, key, value)
		
		self[key] = value
	
	def __delattr__(self, key):
		''' Implement attribute reset using the syntax ``del conf.datapath`` '''

		if key[0] == '_' or key == 'register_variable':
			return default_dict.__delattr__(self, key)
		
		del self[key]
	
	def _add_single_variable(self, q, q_file, q_long, q_units, q_bins):
		''' Add a single variable to the configurarion '''

		if q_file:
			self[self._Q_FILE][q] = q_file
			self[self._Q][q_file] = q
		if q_long:
			self[self._Q_LONG][q] = q_long
		if q_units:
			self[self._Q_UNITS][q] = q_units
		if q_bins:
			self[self._Q_BINS][q] = q_bins
	
	def _unpack_q(self, q_item):
		''' Unpack a tuple holding pertinent information about a variable '''

		if len(q_item) == 2:
			q, q_file = q_item
			qlong = None; q_units = None; q_bins = None
		elif len(q_item) == 3:
			q, q_file, q_long = q_item
			q_units = None; q_bins = None
		elif len(q_item) == 4:
			q, q_file, q_long, q_units = q_item
			q_bins = None
		elif len(q_item) == 5:
			q, q_file, q_long, q_units, q_bins = q_item
		else:
			raise ValueError, 'q_item must be length 2--5, got %d instead.' % len(q_item)

		return q, q_file, q_long, q_units, q_bins


	def register_variable(self, qs, plevs):
		''' Register (a) new variable(s) 
		
		Parameters
		----------
		qs : (list of) 2-,3-,4- or 5-tuple(s)
		    Tuple(s) describing the variable, including the entries 

		        1. Variable name as it appears in the netCDF file
			2. File name segment for files containing the variable
			3. (Optional) Long name for the variable
			4. (Optional) Units of the variale
			5. (Optional) Binning information for the variable

		plevs : list
		    Vertical levels on which the given variables are available. Might be an empty 
		    list if the variable is only to be registered in the non-plot configuration.
		'''

		if type(qs) == list:
			for q in qs:
				q, q_file, q_long, q_units, q_bins = self._unpack_q(q)
				self._add_single_variable(q, q_file, q_long, q_units, q_bins)
			qs = set([q[0] for q in qs])
		else:
			q, q_file, q_long, q_units, q_bins = self._unpack_q(qs)
			self._add_single_variable(q, q_file, q_long, q_units, q_bins)
			qs = set([qs[0],])
		
		if len(plevs) > 0:
			plevs = set(plevs)
			self[self._PLOT].add_table([plevs, qs])
			self[self._PLOTF].add_table([plevs, qs])



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
	'q': {}, 		# Given a file name segment, which variable to we expect to find there?
	'qf': {}, 		# Given a variable name, in which file do we find it?
	'q_units': {},
	'q_long': {},
	'q_bins': {},
	'datapath': ['.', ],
	'opath': '.',
	'plotpath': '.',
	'file_std': '',
	'file_stat': '', 	# Are these general enough to be included? In this case: Would need routines to write those files.
	'file_mstat': '', 	# -"-
	'file_static': None,
	'years': [],
	'times': [],
	'plevs': [],
	'ptlevs': [],
	'pvlevs': [],
	'zlevs': [],
	'mlevs': [],
	'sfclevs': [],
	'plot': plot_settings_dict(PLOT_DEFAULTS),
	'plotf': plot_settings_dict(PLOTF_DEFAULTS),
})

# that's it
