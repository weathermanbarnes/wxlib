#!/usr/bin/env python
# -*- encoding: utf-8


''' Defines default settings and convenient ways to change them

Dynlib settings are in essence a set of key-value pairs, which could be well
represented by the python built-in ``dict`` object. However, the pure ``dict`` object
has some drawbacks in this context:

1. ``dict`` objects cannot represent default values for keys. Whenever a key in
   ``dict`` is overwritten, the original content is lost.
2. ``dict`` cannot represent interdependencies between different keys.
3. ``dict`` cannot represent several variants of itself. As an example, the 
   plot configuration for different variables will be largely identical, but only
   differ in a few configuration keys. When represented by pure ``dict`` objects, 
   the plot configuration for each variable would need to be full independent 
   from each other. This indepence complicates would makes it tedious and 
   error-prone to apply a customised configuration to all variables.

For these reasons, this module introduces a settings_dict object, which is 
subsequently used to define first the default plot configuration for all plots, 
and second some adapted configurarion for specific variables.

Furthermore, this module defines a settings object which contains configuraton
keys as attributes, including the plot configuration objects.

To Do: Potential changes
------------------------

The settings object might be converted into an instance of settings_dict in the 
future. This would allow to use the advantages 1.-3. also for non-plot
configuration.
'''


from copy import copy, deepcopy
import os
import math
from collections import MutableMapping as mutmap

import numpy as np
import matplotlib.pyplot as plt

import proj
import cm




# #############################################################################
# 1. Scales, ticks and labels
# 
scale_oro_c = range(10000,80001,10000)
scale_oro_cf = range(-19000,51000,2000)

scale_pv = np.array([-2,-1, 1, 2])

scale_defang = (np.arange(-18,19)-0.5)*np.pi/36.0
ticks_defang = np.arange(-4,5)*3.1415926535/8.0 
labels_defang = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']

scale_dd = np.arange(0,36.1)*10.0
ticks_dd = np.arange(0,8)*360.0/8.0 
labels_dd = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']

scale_q = np.arange(0.0, 10.1, 1.0)

# #############################################################################
# 2. Default hooks for plotting
#
# TODO: These hooks should be taken into account for both contour and contourf plots automatically, 
#       without the below repetition
hooks = {}
hooks['defabs'] = lambda defabs: defabs*1e5 	# To get typical magnitudes of deformation to the order 1-10
hooks['pv'] = lambda pv: pv*1e6 		# From SI units to PVU
hooks['q'] = lambda q: q*1e3 			# From kg/kg to g/kg
hooks['msl'] = lambda msl: msl/100.0 		# From Pa to hPa
hooks['z'] = lambda msl: msl/98.1		# From m2s-2 to gpdm
def _tmp(oro):
	oro[oro <= 100] = -17000
	return oro
hooks['oro'] = _tmp



# #############################################################################
# 3. Default settings
#

Q = {'defabs': 'defabs', 'defang': 'defang', 'm': 'mont', 'p': 'pres', 'msl': 'msl', 'u': 'u', 'v': 'v', 'w': 'w', 'q': 'q', 
		'T': 't', 'th': 'pt', 'the': 'thetae', 'Z': 'z', 'oro': 'oro', 'rsr': 'rsr', 'ow': 'ow', 'pv': 'pv', 
		'front': 'front', 'cold_front': 'cold_front', 'warm_front': 'warm_front', 'stat_front': 'stat_front',
		'convl': 'convl', 'defl': 'defl', 'vorl': 'vorl', 
		'jetaxis': 'jetaxis', 'tw': 'tcw', 'wv': 'tcwv', 'zeta': 'vo', 'div': 'div', 'ps': 'sp', 'ff': 'ff',
		'ttr': 'ttr', 'SST': 'sst'}
QI = {'defabs': 'defabs', 'defang': 'defang', 'mont': 'm', 'pres': 'p', 'msl': 'msl', 'u': 'u', 'v': 'v', 'w': 'w', 'q': 'q', 
		't': 'T', 'pt': 'th', 'thetae': 'the', 'z': 'Z', 'oro': 'oro', 'rsr': 'rsr', 'ow': 'ow', 'pv': 'pv', 
		'front': 'front', 'cold_front': 'cold_front', 'warm_front': 'warm_front', 'stat_front': 'stat_front',
		'froff': 'front', 'convl': 'convl', 'cloff': 'convl', 
		'defl': 'defl', 'dloff': 'defl', 'vorl': 'vorl', 'vloff': 'vorl', 
		'jetaxis': 'jetaxis', 'jaoff': 'jetaxis', 'tcw': 'tw', 'tcwv': 'wv', 'vo': 'zeta',
		'div': 'div', 'sp': 'ps', 'ff': 'ff',
		'ttr': 'ttr', 'sst': 'SST'}

UNITS = {'defabs': 's-1', 'defang': 'rad', 'defanr': 'rad', 'the': 'K', 'rsr': '1', 'ow': 's-2', 
		'zeta': 's-1', 'div': 's-1', 'ttr': 'W m-2' }
LONG = {'defabs': 'Total deformation', 'defang': 'Deformation angle', 'defanr': 'Deformation angle in natural coordinates',
		'the': 'Equivalent potential temperature', 'rsr': 'Rotation/Strain-ratio', 'ow': 'Okubo-Weiss criterion',
		'zeta': 'Horizontal vorticity', 'div': 'Horizontal divergence', 'jetaxis': 'Jet axis lines',
		'ff': 'wind speed', 'pv': 'Potential vorticity', 'msl': 'Mean sea level pressure', 'w': 'Vertical velocity',
		'cold_front': 'Cold front lines', 'warm_front': 'Warm front lines', 'stat_front': 'Stationary front lines',
		'ttr': 'Top net thermal radiation',
}

_rose = [17,]
_rose.extend(range(-18,18))
BINS_Q = {'defang': np.array(_rose)*math.pi/36.0+math.pi/72.0, }

DATAPATH = ['.', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY', '/Data/gfi/users/csp001/share', ]
OPATH    = '.'
PPATH    = '.'
FILE_STD   = 'ei.ans.%(time)d.%(plev)s.%(q)s'
FILE_STAT  = 'ei.ans.%(time)d.%(plev)s.%(q)s.stat'
FILE_MSTAT = 'ei.ans.stat.%(plev)s.%(q)s'
STD_SLICE  = (slice(None), slice(None), slice(None))
YEARS  = range(1979,2015)
PLEVS  = ['100', '200', '300', '400', '500', '550', '600', '650', '700', '750', '800', '850', '900', '950', '1000', ]
PTLEVS = ['pt300', 'pt315', 'pt330', 'pt350', ]
PVLEVS = ['pv2000', ]

# DEFAULT contour settings
if os.getenv('DYNLIB_PLOT_PRINT'):
	DEFAULT_KWARGS = {'m': proj.world, 'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 'auto', 
		'overlays': [], 'disable_cb': True, 'show': False, 'save': '', 'title': '', 'hook': None,
		'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroscale': scale_oro_c,
		'oroalpha': 0.4, 'ticks': None, 'ticklabels': [], 'scale_symmetric_zero': False}
else:
	DEFAULT_KWARGS = {'m': proj.world, 'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 'auto', 
		'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 'hook': None,
		'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroscale': scale_oro_c,
		'oroalpha': 0.4, 'ticks': None, 'ticklabels': [], 'scale_symmetric_zero': False }
DEFAULT_CONTOUR_KWARGS = {'colors': 'k', 'alpha': 1.0, 'cmap': None, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'neither', 'linewidths': 2.0, 'linestyles': None, 'contour_labels': False}
DEFAULT_CONTOURF_KWARGS = {'colors': None, 'alpha': 1.0, 'cmap': cm.gist_ncar, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'both', 'hatches': None }

MUTEX_GROUPS = [set(['colors', 'cmap']), ]

# DEFAULT settings per quantity Q on contourf plots
DEFAULT_Q_C = {}
DEFAULT_Q_C['defabs'] = {'hook': hooks['defabs']}
DEFAULT_Q_C['pv']  = {'hook': hooks['pv']}
DEFAULT_Q_C['q'] = {'hook': hooks['q']}
DEFAULT_Q_C['msl'] = {'hook': hooks['msl']}
DEFAULT_Q_C['z'] = {'hook': hooks['z']}

# DEFAULT settings per quantity Q on contourf plots
DEFAULT_Q_CF = {}
DEFAULT_Q_CF['defabs'] = {'cmap': cm.defabs(), 'hook': hooks['defabs']}
DEFAULT_Q_CF['defang'] = {'cmap': cm.periodic(), 'scale': scale_defang, 'ticks': ticks_defang, 
	'ticklabels': labels_defang}
DEFAULT_Q_CF['t']   = {'cmap': plt.cm.RdBu_r}
DEFAULT_Q_CF['pt']   = {'cmap': plt.cm.RdBu_r}
DEFAULT_Q_CF['thetae']   = {'cmap': plt.cm.RdBu_r}
DEFAULT_Q_CF['pv']  = {'hook': hooks['pv']}
DEFAULT_Q_CF['oro'] = {'scale': scale_oro_cf, 'cmap': plt.cm.gist_earth, 'hook': hooks['oro']}
DEFAULT_Q_CF['q'] = {'cmap': cm.q(), 'hook': hooks['q']}
DEFAULT_Q_CF['tcw'] = {'cmap': cm.q()}
DEFAULT_Q_CF['tcwv'] = {'cmap': cm.q()}
DEFAULT_Q_CF['msl'] = {'hook': hooks['msl']}
DEFAULT_Q_CF['z'] = {'hook': hooks['z']}


# #############################################################################
# 4. Making the settings easily available
#
class settings_dict(mutmap):
	''' An extended dictionary object to be used for configuration
	
	It handles default values:

	 * They are readable through the ``default`` attribute.
	 * They can be restored by using the ``reset`` function.
	
	It handles interdependencies between keys:
	
	 * Mutually exclusive keys can be defined by filling the
	   ``_mutex_group`` attribute during the initialisation.
	
	Furthermore, it provides convenience functions for working with the 
	configuration and defines the neccessary magic functions to appear 
	largely like a pure python ``dict``.
	'''
	_mutexes = {}

	def __init__(self, init={}, contourobj=None):
		''' Initialisation
		
		Parameters
		----------
		init : dict
			Optional initial overrides to the default configuration
		contourobj : settings_dict
			Optional. Take default values from this configuration object. 
		'''
		for mutex_group in MUTEX_GROUPS:
			for key in mutex_group:
				self._mutexes[key] = copy(mutex_group)
				self._mutexes[key].remove(key)

		if contourobj:
			self.default = contourobj.default
		else: 
			self.default = None
		self.default_q = copy(init)
		self._ = init

		mutmap.__init__(self)

		return

	
	def __getitem__(self, key):
		''' Implements access through ``conf['key']`` '''
		if key not in self._ and self.default:
			return self.default[key]
		else: 	
			return self._[key]
	

	def __setitem__(self, key, value):
		''' Implements setting overrides through ``conf['key'] = value`` '''
		if self.default and key not in self.default:
			raise KeyError, 'Cannot add new keys to the configuration'
		#elif key not in self:
		#	raise KeyError, 'Cannot add new keys to the configuration'

		if value != None:
			for mutex in self._mutexes.get(key, []):
				self._[mutex] = None
		self._[key] = value

		return


	def __delitem__(self, key):
		''' Implements deleting overrides by ``del conf['key']`` '''
		self.reset(key)

		return
	

	def __contains__(self, key):
		''' Implements checks with the syntax ``if key in conf`` '''
		return key in self.keys()


	def __iter__(self):
		''' Implements loops over the configuration keys, syntax ``for key in conf`` '''
		if self.default:
			for key in self.default:
				yield key
		else:
			for key in self._:
				yield key
	
	iterkeys = __iter__

	#def __nonzero__(self):
	#	return True


	def __eq__(self, other):
		''' Check if two objects are equal '''
        	return sorted(self.iteritems()) == sorted(other.iteritems())



	def __len__(self):
		''' Return a length of the object, here number of defined configuration keys '''
		if self.default:
			return len(self.default)
		else:
			return len(self._)
	

	def __repr__(self):
		''' Give a string-representation of the object '''
		return dict(self.iteritems()).__repr__()

	
	def __reset_key(self, key):
		''' Reset a single key to its defaults '''
		if key in self.default_q:
			self._[key] = self.default_q[key]
		else:
			del self._[key]

		return


	def iteritems(self):
		''' Implements loops over the configuration keys, syntax ``for key, value in conf.iteritems()`` '''
		for key in self.iterkeys():
			yield key, self[key]
	

	#def items(self):
	#
	#	return


	def values(self):
		''' Stub method to implements loops over the configuration keys, syntax ``for value in conf.values()`` 
		
		This method is not implemented, but only present for consistency 
		with the API of pure python ``dict`` objects.
		'''
		raise NotImplementedError, 'If you need it, implement it!'
	itervalues = values
	viewitems = values
	viewvalues = values
	

	#def get(self, key):
	#	if key in self._:
	#		return self._[key]
	#	else:
	#		return default
	
	def merge(self, **kwargs):
		''' Merge the given keyword arguments with the configuration in this object

		Keyword arguments
		-----------------
		Anything

		Returns
		-------
		dict
		    A dictionary containing the configuration in this object, potentially overriden
		    and extended by the given keyword arguments.
		'''

		rkwargs = dict(self)
		for kwarg, argv in kwargs.items():
			rkwargs[kwarg] = argv

		return rkwargs

	def reset(self, key=None):
		''' Reset the configuration to its defaults (optional for a specific key)

		Parameters
		----------
		key : any valid dict key
		    Optional: Key to be reset to the defaults.
		'''
		if not self.default:
			raise TypeError, 'Reset not possible for default object!'

		if key:
			for mutex in self._mutexes.get(key, []):
				if mutex in self._:
					self.__reset_key(mutex)
			self.__reset_key(key)
		else:
			self._ = copy(self.default_q)

		return



class settings_contour(object):
	''' Container object for contour plot settings

	The object contains pointers to the respective settings for each
	registered variable, available through attributes of this object.
	''' 

	default = settings_dict(DEFAULT_CONTOUR_KWARGS)
	default.update(DEFAULT_KWARGS)
	default_q = DEFAULT_Q_C

	_overrides = {}

	def __init__(self):
		for q in QI:
			self.new(q, self.default_q.get(q, {}))

		return


	def __getattribute__(self, q):
		''' Implement configuration access through the syntax ``conf.key`` '''

		if q[0] == '_':
			return object.__getattribute__(self, q)
		elif q not in self._overrides:
			return object.__getattribute__(self, q)
		
		return self._overrides[q]


	def __setattr__(self, q, value):
		''' Implement configuration cjamhes through the syntax ``conf.key = value`` '''

		if q in self._overrides or q == 'default':
			raise AttributeError, 'The attributes cannot be overwritten'
		object.__setattr__(self, q, value)

		return


	def new(self, q, conf):
		''' Register new variable 

		Parameters
		----------
		q : str
		    Abbreviated variable name, as it would appear in as a netCDF variable name
		conf : dict-like
		    Initial overrides for the plot configuration of the new variable
		'''
		self._overrides[q] = settings_dict(conf, self)

		return
	

	def merge(self, q, **kwargs):
		''' Proxy for ``settings_dict.merge`` 
		
		Parameters
		----------
		q : str
		    Abbreviated variable name, as it would appear in as a netCDF variable name

		Keyword Arguments
		-----------------
		Anything
		'''

		return self.__getattribute__(q).merge(**kwargs)


	def reset(self, q, key=None):
		''' Reset the configuration to its defaults (optional for a specific key)

		Parameters
		----------
		key : any valid dict key
		    Optional: Key to be reset to the defaults.
		'''
		self.__getattribute__(q).reset(key)

		return



class settings_contourf(settings_contour):
	''' Container object for filled contour plot settings

	The object contains pointers to the respective settings for each
	registered variable, available through attributes of this object.
	''' 

	default = settings_dict(DEFAULT_CONTOURF_KWARGS)
	default.update(DEFAULT_KWARGS)
	default_q = DEFAULT_Q_CF
	_overrides = {}
	


# TODO: Should this derive from settings_dict? Or should it be a settings_dict instance?
# - Would prohibit (accientally but also volutarily) setting new config keys
# - Would potentially avoid some code duplication: __special_funcs__() and reset()

class settings(object):
	''' Container object for all settings

	Note: This object might become a settings_dict instance in the future!
	
	The object is yet another interface to a dictionary, including default values.
	In contrast to the ``settings_dict`` object, here the configuration keys are
	made accessible through attributes, with the syntax being ``conf.key`` 
	instead of the dictionary lookup syntax ``conf['key']``.
	'''

	__default = {
		'q': Q,
		'qi': QI,
		'q_units': UNITS,
		'q_long': LONG,
		'bins': BINS_Q,
		'datapath': DATAPATH,
		'opath': OPATH,
		'ppath': PPATH,
		'file_std': FILE_STD,
		'file_stat': FILE_STAT,
		'file_mstat': FILE_MSTAT,
		'std_slice': STD_SLICE,
		'years': YEARS,
		'plevs': PLEVS,
		'ptlevs': PTLEVS,
		'pvlevs': PVLEVS,
		'contour': settings_contour(),
		'contourf': settings_contourf(),
	}

	def __init__(self):
		self.__current = deepcopy(self.__default)

		return


	def __getattribute__(self, key):
		''' Implement configuration access through the syntax ``conf.key`` '''

		if key[0] == '_' or key in ['reset', 'new_variable', ]:
			return object.__getattribute__(self, key)

		return self.__current[key]


	def __setattr__(self, key, value):
		''' Implement configuration cjamhes through the syntax ``conf.key = value`` '''

		if key in ['_settings__default', '__setattr__']:
			raise AttributeError, 'The default values cannot be overwritten'
		elif key[0] == '_' or key in ['reset', ]:
			return object.__setattr__(self, key, value)
		self.__current[key] = value

		return
	

	def new_variable(self, q, qlong, conf_cf={}, conf_c={}):
		''' Register a new variable in the plot and metadata configuration

		Parameters
		----------
		q : str
		    Abbreviated variable name, as it would appear in as a netCDF variable name
		qlong: str
		    Full variable name, as it would appear in the long_name attribute of netCDF variables
		conf_cf : dict-like
		    Optional. Standard configuration overrides for filled-contour plots of the new variable
		conf_c : dict-like
		    Optional. Standard configuration overrides for contour plots of the new variable
		'''
		self.contourf.new(qlong, conf_cf)
		self.contour.new(qlong, conf_c)
		self.q[q] = qlong
		self.qi[qlong] = q

		return


	def reset(self, key=None):
		''' Reset the configuration to its defaults (optional for a specific key)

		Parameters
		----------
		key : any valid dict key
		    Optional: Key to be reset to the defaults.
		'''
		if key and key in self.__current:
			self.__current[key] = copy(self.__default[key])
		else:
			self.__current = deepcopy(self.__default)

		return

conf = settings()



# #############################################################################
# 5. Clean-Up: Making the default settings only available through settings objects
# 
del Q, QI, BINS_Q, UNITS, LONG, DATAPATH, OPATH, PPATH, FILE_STD, FILE_STAT, FILE_MSTAT, STD_SLICE, YEARS, PLEVS, PTLEVS, PVLEVS
del DEFAULT_KWARGS, DEFAULT_CONTOUR_KWARGS, DEFAULT_CONTOURF_KWARGS, DEFAULT_Q_C, DEFAULT_Q_CF


# that's it
