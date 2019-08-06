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
 #. ``settings_obj``: A variant of the ``default_dict`` where access is also possible
    through the attribute syntax dat.key. This class is the basis for the conf 
    object, the configuration root. Within, ``conf.plot`` and ``conf.plotf`` are 
    ``plot_settings_dict`` objects.
 #. ``in_context``, ``def_context`` and ``get_context`` functions, to define new 
    contexts (i.e. named modular sets of settings), to check whether a specific 
    context is active, and to access the settings associated with a context, 
    respectively.
'''

from __future__ import absolute_import, unicode_literals

import numpy as np
from copy import copy 

import collections


__context__ = dict()

def in_context(*contexts):
	''' Check if any of the given contexts is loaded
	
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

def def_context(context, parent=None):
	''' Define a new context

	The new context will contain a plot configuration if 
	 - the default config does contain one when creating a new context
	 - the parent context contains one when extending another context

	Parameters
	----------
	context : str
	    Name of the context to be defined
	parent : str
	    Optional: Name of the parent context or the string "active". If given, no 
	    new set of configuration is created, and the context is expected to operate 
	    on the same configuration as its parent. If "active", whatever context is 
	    currently active is used as parent. If no parent given, the newly created 
	    context is made active immediately.
	
	Returns
	-------
	settings_obj
	    The newly created or parent context configuration
	'''

	if parent == 'active':
		conf_ = conf
	elif parent:
		conf_ = __context__[parent]
	else:
		_plot = {}
		_plotf = {}
		if not context == 'default':
			default_ = get_context('default')
			if not type(default_.plot) == dict:
				from .settings import plot_settings_dict, PLOT_DEFAULTS, PLOTF_DEFAULTS, MUTEX_GROUPS

				_plot = plot_settings_dict(PLOT_DEFAULTS, MUTEX_GROUPS)
				_plotf = plot_settings_dict(PLOTF_DEFAULTS, MUTEX_GROUPS)

		conf_ = settings_obj({
			'q': {}, 		# Given a file name segment, which variable to we expect to find there?
			'qf': {}, 		# Given a variable name, in which file do we find it?
			'q_units': {},
			'q_long': {},
			'q_bins': {},
			'datapath': ['.', '/Data/gfi/users/local/share'], 
			'opath': '.',
			'oformat': 'nc',
			'plotpath': '.',
			'plotformat': 'png',
			'file_std': None,
			'file_agg': None,
			'file_ts': None,
			'file_timeless': None,
			'file_static': None,
			'epoch': None,
			'years': [],
			'times': [],
			'timestep': None,
			'gridsize': None,
			'local_timezone': 'Europe/Oslo',
			'plevs': [],
			'ptlevs': [],
			'pvlevs': [],
			'zlevs': [],
			'mlevs': [],
			'sfclevs': [],
			'plot': _plot, 
			'plotf': _plotf, 
		}, [])
		
	__context__[context] = conf_

	if not parent:
		set_active_context(context)
	
	return conf_

def get_context(context):
	''' Get the configuration for the given context
	
	Parameters
	----------
	context : str
	    Context queried for configuration
	
	Returns
	-------
	settings_obj
	    The currently active configuration
	'''
	
	return __context__[context]

def set_active_context(context):
	''' Set the currently active configuration by context name 
	
	Parameters
	----------
	context : str
	    Context to be made active
	'''
	
	global conf
	conf = __context__[context]

def is_active_context(context):
	''' Query if the given context is currently active
	
	Parameters
	----------
	context : str
	    Context to be checked
	
	Returns
	-------
	bool
	    ``True'' if the given context is active
	'''
	
	global conf
	return id(conf) == id(__context__[context])

def get_active_context():
	''' Get the currently active configuration
	
	Returns
	-------
	settings_obj
	    The currently active configuration
	'''
	
	global conf
	return conf


class default_dict(collections.MutableMapping):
	''' Dictionary with a predefined set of valid keys and default values for each 

	Non-default values are stored in ``self._``, the defaults in ``self._defaults``.
	
	Parameters
	----------
	defaults : dict
	    Definition of valid keys and their default values.
	mutexes : list of sets
	    Sets of mutually exclusive configuration keys. As soon as one is not None, all others 
	    are automaticall set to None.
	'''

	def __init__(self, defaults, mutexes):
		self._ = {}
		self._defaults = defaults

		self._mutexes = {}
		for mutex in mutexes:
			for key in mutex:
				self._mutexes[key] = copy(mutex)
				self._mutexes[key].remove(key)

		# Mutable (approximated by "non-hashable") types can be changed in place. 
		# Hence, they must be copied over to the overrides dict to avoid overriding defaults
		for key, value in defaults.items():
			if not isinstance(value, collections.Hashable):
				self._[key] = copy(value)

		collections.MutableMapping.__init__(self)
		
	def __getitem__(self, key):
		''' Implements the access syntax ``dat[key]`` '''

		if key in self._:
			return self._[key]
		elif key in self._defaults:
			return self._defaults[key]
		else:
			raise KeyError(key)

	def __setitem__(self, key, value):
		''' Implements the assignment syntax ``dat[key] = value`` '''
		
		# Respect that some combination configuration are not meaningful
		if not type(value) == type(None):
			for mutex in self._mutexes.get(key, []):
				self._[mutex] = None

		if key in self._defaults:
			self._[key] = value
		else:
			raise KeyError(key)

	def __delitem__(self, key):
		''' Implements the reset syntax ``del dat[key]`` '''

		if key in self._:
			# Mutable (approximated by "non-hashable") types can be changed in place. 
			# Hence, they must be copied over to the overrides dict to avoid overriding defaults
			if not isinstance(self._defaults[key], collections.Hashable):
				self._[key] = copy(self._defaults[key])
			else:
				del self._[key]
		elif key in self._defaults:
			pass
		else:
			raise KeyError(key)

	def __iter__(self):
		''' Implement interation over the dictionary, e.g. the syntax ``for key in dat:`` '''

		return self._defaults.__iter__()

	def __len__(self):
		''' Implements the length query, allowing among others to use ``len(dat)`` '''

		return len(self._defaults)


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
	    Definition of valid keys and their default values.
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
			raise TypeError('_get() requires at least one key')
	
	def _add_default(self, key, value):
		''' Add a new configuration item and its default value '''

		if key in self._defaults:
			raise KeyError(str(key) + ' already set')
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

		if not type(q_file) == type(None):
			self[self._Q_FILE][q] = q_file
			self[self._Q][q_file] = q
		if q_long:
			self[self._Q_LONG][q] = q_long
		if q_units:
			self[self._Q_UNITS][q] = q_units
		if type(q_bins) in [list, np.ndarray]:
			self[self._Q_BINS][q] = q_bins
	
	def _unpack_q(self, q_item):
		''' Unpack a tuple holding pertinent information about a variable '''

		if len(q_item) == 2:
			q, q_file = q_item
			q_long = None; q_units = None; q_bins = None
		elif len(q_item) == 3:
			q, q_file, q_long = q_item
			q_units = None; q_bins = None
		elif len(q_item) == 4:
			q, q_file, q_long, q_units = q_item
			q_bins = None
		elif len(q_item) == 5:
			q, q_file, q_long, q_units, q_bins = q_item
		else:
			raise ValueError('q_item must be length 2--5, got %d instead.' % len(q_item))

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
		    list if the variable is only to be registered for a general level "None".
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
		
		plevs = set(plevs)
		
		# Optionally: Add plot configuration for the new variable
		if hasattr(self[self._PLOT], 'add_table'):
			self[self._PLOT].add_table([copy(plevs), copy(qs)])
		if hasattr(self[self._PLOTF], 'add_table'):
			self[self._PLOTF].add_table([copy(plevs), copy(qs)])


conf = def_context('default')

# C'est la fin.
