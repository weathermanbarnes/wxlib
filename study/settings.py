#!/usr/bin/env python
# -*- encoding: utf-8


''' Defines basic data structures for settings storage and convenient access

To be coded and written.
'''

from collections import MutableMapping as mutmap


__context__ = set([])

def in_context(*contexts):
	for c in contexts:
		if c in __context__:
			return True
	
	return False

def def_context(context):
	__context__.add(context)



def register_variable(q, q_file=None, q_long='', q_unit=''):
	pass

def register_level(plev):
	pass


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
	
	The set of allowed keys for each dimension is limited, and mainly set during 
	initialisation. Valid keys can only be introduced later through the ``add_default()``
	and ``add_in_dimension()`` functions.
	"""
	def __init__(self, first_dimensions, defaults):
		self._fdims = first_dimensions 		# List of sets containing valid keys for each dimension
		self._ndims = len(self._fdims) + 1
		default_dict.__init__(self, defaults) 	# dictionary of valid keys in the last dimension and their default values
	
	def __getitem__(self, key):
		if len(key) < self._ndims - 1:
			raise KeyError, 'Key must contain at least %d dimensions' % len(self._fdims)
		elif len(key) > self._ndims:
			raise KeyError, 'Key can at most contain %d dimensions' % (len(self._fdims) + 1)
		
		for i, valid_skeys in zip(range(len(self._fdims)), self._fdims):
			if not key[i] in valid_skeys:
				raise KeyError, key
		
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


class hierarchical_default_dict(default_dict):
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

# that's it
