#!/usr/bin/env python
# -*- encoding: utf-8

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
hooks = {}
hooks['defabs'] = lambda defabs: defabs*1e5
hooks['pv']     = lambda pv: pv*1e6
hooks['q']     = lambda q: q*1e3
def _tmp(oro):
	oro[oro <= 100] = -17000
	return oro
hooks['oro'] = _tmp



# #############################################################################
# 3. Default settings
#
Q = {'defabs': 'defabs', 'defang': 'defang', 'm': 'mont', 'p': 'pres', 'u': 'u', 'v': 'v', 'q': 'q', 
		'T': 't', 'the': 'thetae', 'Z': 'z', 'oro': 'oro', 'rsr': 'rsr', 'ow': 'ow', 'pv': 'pv', 
		'fronts': 'fronts', 'convls': 'convls', 'defls': 'defls', 'vorls': 'vorls', 
		'jetaxis': 'jetaxis', 'tw': 'tcw', 'wv': 'tcwv', 'zeta': 'vo', 'div': 'div', 'ps': 'sp'}
QI = {'defabs': 'defabs', 'defang': 'defang', 'mont': 'm', 'pres': 'p', 'u': 'u', 'v': 'v', 'q': 'q', 
		't': 'T', 'thetae': 'the', 'z': 'Z', 'oro': 'oro', 'rsr': 'rsr', 'ow': 'ow', 'pv': 'pv', 
		'fronts': 'fronts', 'froff': 'fronts', 'convls': 'convls', 'cloff': 'convls', 
		'defls': 'defls', 'dloff': 'defls', 'vorls': 'vorls', 'vloff': 'vorls', 
		'jetaxis': 'jetaxis', 'jaoff': 'jetaxis', 'tcw': 'tw', 'tcwv': 'wv', 'vo': 'zeta',
		'div': 'div', 'sp': 'ps'}

UNITS = {'defabs': 's-1', 'defang': 'rad', 'defanr': 'rad', 'the': 'K', 'rsr': '1', 'ow': 's-2', 
	'zeta': 's-1', 'div': 's-1'}
LONG = {'defabs': 'Total deformation', 'defang': 'Deformation angle', 'defanr': 'Deformation angle in natural coordinates',
		'the': 'Equivalent potential temperature', 'rsr': 'Rotation/Strain-ratio', 'ow': 'Okubo-Weiss criterion',
		'zeta': 'Horizontal vorticity', 'div': 'Horizontal divergence', 'jetaxis': 'Jet axis lines'}

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
YEARS  = range(1979,2013)
PLEVS  = ['100', '200', '300', '400', '500', '550', '600', '650', '700', '750', '800', '850', '900', '950', '1000', ]
PTLEVS = ['pt300', 'pt315', 'pt330', 'pt350', ]
PVLEVS = ['pv2000', ]

# DEFAULT contour settings
if os.getenv('DYNLIB_PLOT_PRINT'):
	DEFAULT_KWARGS = {'m': proj.wmap, 'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 'auto', 
		'overlays': [], 'disable_cb': True, 'show': False, 'save': '', 'title': '', 'hook': None,
		'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroscale': scale_oro_c,
		'oroalpha': 0.4, 'ticks': None, 'ticklabels': [], 'scale_symmetric_zero': False}
else:
	DEFAULT_KWARGS = {'m': proj.wmap, 'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 'auto', 
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

# DEFAULT settings per quantity Q on contourf plots
DEFAULT_Q_CF = {}
DEFAULT_Q_CF['defabs'] = {'cmap': cm.defabs2(), 'hook': hooks['defabs']}
DEFAULT_Q_CF['defang'] = {'cmap': cm.periodic3(), 'scale': scale_defang, 'ticks': ticks_defang, 
	'ticklabels': labels_defang}
DEFAULT_Q_CF['t']   = {'cmap': plt.cm.RdBu_r}
DEFAULT_Q_CF['thetae']   = {'cmap': plt.cm.RdBu_r}
DEFAULT_Q_CF['pv']  = {'hook': hooks['pv']}
DEFAULT_Q_CF['oro'] = {'scale': scale_oro_cf, 'cmap': plt.cm.gist_earth, 'hook': hooks['oro']}
DEFAULT_Q_CF['q'] = {'cmap': cm.q(), 'hook': hooks['q']}
DEFAULT_Q_CF['tcw'] = {'cmap': cm.q()}
DEFAULT_Q_CF['tcwv'] = {'cmap': cm.q()}


# #############################################################################
# 4. Making the settings easily available
#
class settings_dict(mutmap):
	_mutexes = {}

	def __init__(self, init={}, contourobj=None):
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
		if key not in self._ and self.default:
			return self.default[key]
		else: 	
			return self._[key]
	

	def __setitem__(self, key, value):
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
		self.reset(key)

		return
	

	def __contains__(self, key):
		return key in self.keys()


	def __iter__(self):
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
        	return sorted(self.iteritems()) == sorted(other.iteritems())



	def __len__(self):
		if self.default:
			return len(self.default)
		else:
			return len(self._)
	

	def __repr__(self):
		return dict(self.iteritems()).__repr__()

	
	def __reset_key(self, key):
		if key in self.default_q:
			self._[key] = self.default_q[key]
		else:
			del self._[key]

		return


	def iteritems(self):
		for key in self.iterkeys():
			yield key, self[key]
	

	#def items(self):
	#
	#	return


	def values(self):
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
		rkwargs = dict(self)
		for kwarg, argv in kwargs.items():
			rkwargs[kwarg] = argv

		return rkwargs

	def reset(self, key=None):
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
	default = settings_dict(DEFAULT_CONTOUR_KWARGS)
	default.update(DEFAULT_KWARGS)
	default_q = DEFAULT_Q_C

	_overrides = {}

	def __init__(self):
		for q in QI:
			self.new(q, self.default_q.get(q, {}))

		return


	def __getattribute__(self, q):
		if q[0] == '_':
			return object.__getattribute__(self, q)
		elif q not in self._overrides:
			return object.__getattribute__(self, q)
		
		return self._overrides[q]


	def __setattr__(self, q, value):
		if q in self._overrides or q == 'default':
			raise AttributeError, 'The attributes cannot be overwritten'
		object.__setattr__(self, q, value)

		return


	def new(self, q, conf):
		self._overrides[q] = settings_dict(conf, self)

		return
	

	def merge(self, q, **kwargs):
		return self.__getattribute__(q).merge(**kwargs)


	def reset(self, q, key=None):
		self.__getattribute__(q).reset(key)

		return



class settings_contourf(settings_contour):
	default = settings_dict(DEFAULT_CONTOURF_KWARGS)
	default.update(DEFAULT_KWARGS)
	default_q = DEFAULT_Q_CF
	_overrides = {}
	


class settings(object):
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
		if key[0] == '_' or key in ['reset', 'new_variable', ]:
			return object.__getattribute__(self, key)

		return self.__current[key]


	def __setattr__(self, key, value):
		if key in ['_settings__default', '__setattr__']:
			raise AttributeError, 'The default values cannot be overwritten'
		elif key[0] == '_' or key in ['reset', ]:
			return object.__setattr__(self, key, value)
		self.__current[key] = value

		return
	

	def new_variable(self, q, qlong, conf_cf={}, conf_c={}):
		self.contourf.new(qlong, conf_cf)
		self.contour.new(qlong, conf_c)
		self.q[q] = qlong
		self.qi[qlong] = q

		return


	def reset(self, key=None):
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
