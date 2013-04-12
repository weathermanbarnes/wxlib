#!/usr/bin/env python
# -*- encoding: utf-8

from copy import copy, deepcopy
import os
import math
from collections import MutableMapping as mutmap

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from proj import wmap



# #############################################################################
# 1. Colour maps
# 
def _get_greys_cm():
	cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'green': ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'blue':  ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1))  }

	return mpl.colors.LinearSegmentedColormap('my_grey',cdict,256)

def _get_defabs_cm():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.4, 0.4), (0.867, 1.0, 1.0), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.5, 0.5), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 1.0, 1.0), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def _get_defabs_cm2():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def _get_q_cm():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.30, 0.30),  (0.867, 0.1, 0.1), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.65, 0.65),  (0.867, 0.2, 0.2), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 0.80, 0.80),  (0.867, 0.6, 0.6), (1.0, 0.8, 0.8))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def _get_periodic_cm():
	cdict = {'red':   ((0.0, 0.0, 0.0), (0.25, 0.8, 0.8), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0)),
		 'green': ((0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.9, 0.9), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 0.6, 0.6), (0.25, 0.0, 0.0), (0.5, 0.2, 0.2), (0.75, 0.0, 0.0), (1.0, 0.6, 0.6))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic',cdict,256)

def _get_periodic_cm2():
	cdict = {'red':   ((0.0, 0.3, 0.3), (0.25, 0.4, 0.4), (0.5, 0.6, 0.6), (0.75, 0.8, 0.8), (1.0, 0.3, 0.3)),
		 'green': ((0.0, 0.2, 0.2), (0.25, 0.8, 0.8), (0.5, 0.6, 0.6), (0.75, 0.4, 0.4), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 0.5, 0.5), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.5, 0.5))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic2',cdict,256)


def _get_periodic_cm3():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.03, 0.9, 0.9), (0.27, 0.3, 0.2), (0.515, 0.5, 0.5), (0.76, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.03, 0.9, 0.9), (0.27, 0.3, 0.2), (0.515, 0.5, 0.5), (0.76, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.03, 1.0, 1.0), (0.27, 0.7, 0.7), (0.515, 0.5, 0.5), (0.76, 0.20, 0.3), (1.0, 0.9, 1.0)) }

	return mpl.colors.LinearSegmentedColormap('my_periodic3',cdict,256)




# #############################################################################
# 2. Scales, ticks and labels
# 
scale_oro = range(10000,80001,10000)
scale_oro_full = range(-19000,51000,2000)

scale_Z_diff = np.arange(-5250,5251,500)

scale_u     = np.arange(20,71,10)
scale_u_diff = np.arange(-30,31,5)

scale_pv = np.array([-2,-1, 1, 2])

scale_defabs = np.arange(5.0,30.1,5.0)
scale_defabs_mean = np.arange(2.0,12.1,2.0)

scale_defang = (np.arange(-18,19)-0.5)*np.pi/36.0
scale_defang_coarse = np.arange(-4,5)*np.pi/8.0 - np.pi/72.0
ticks_defang = np.arange(-4,5)*3.1415926535/8.0 
labls_defang = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']

scale_ow_mean = np.arange(-0.35, 0.36, 0.1)*1.0e-9

scale_pvstir_mean = np.arange(-2.25, 2.26, 0.5)*1.0e-6
scale_pvfil_mean  = np.arange(-7.0, 7.1, 2.0)*1.0e-6

scale_rsr_mean = np.arange(-13.5, 13.6, 3.0)

scale_q = np.arange(0.0, 10.1, 1.0)*1.0e-3
scale_qfs = np.arange(-1.65,1.66,0.3)*1.0e-4
scale_qfs_mean = np.arange(-4.5,4.6,1.0)*1.0e-5

# #############################################################################
# 3. Default hooks for plotting
#
hooks = {}
hooks['defabs'] = lambda defabs: defabs*1e5
hooks['pv']     = lambda pv: pv*1e6
def _tmp(oro):
	oro[oro <= 100] = -17000
	return oro
hooks['oro'] = _tmp



# #############################################################################
# 4. Default settings
#
Q = {'defabs': 'defabs', 'defang': 'defang', 'm': 'mont', 'p': 'pres', 'u': 'u', 'v': 'v', 'q': 'q', 'qstir': 'qstir', 'qfil': 'qfil',
		'T': 't', 'the': 'thetae', 'thefil': 'thetaefil', 'thestir': 'thetaestir', 
		'Z': 'z', 'oro': 'oro', 'rsr': 'rsr', 'ow': 'ow', 'pv': 'pv', 'pvstir': 'pvstir', 'pvfil': 'pvfil', }
_rose = [17,]
_rose.extend(range(-18,18))
BINS_Q = {'defang': np.array(_rose)*math.pi/36.0+math.pi/72.0, }

DATAPATH = ['.', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY', '/Data/gfi/users/csp001/share', ]
OPATH    = '.'
FILE_STD   = 'ei.ans.%d.%s.%s'
FILE_STAT  = 'ei.ans.%d.%s.%s.stat'
FILE_MSTAT = 'ei.ans.stat.%s.%s'
STD_SLICE  = (slice(None), slice(None), slice(None))
YEARS  = range(1979,2012)
PLEVS  = ['100', '200', '300', '400', '500', '550', '600', '650', '700', '750', '800', '850', '900', '950', '1000', ]
PTLEVS = ['pt300', 'pt315', 'pt330', 'pt350', ]
PVLEVS = ['pv2000', ]

# DEFAULT contour settings
if os.getenv('DYNLIB_PLOT_PRINT'):
	DEFAULT_KWARGS = {'m': wmap, 'plev': '800', 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
		'overlays': [], 'disable_cb': True, 'show': False, 'save': '', 'title': '', 'hook': None,
		'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroscale': scale_oro,
		'oroalpha': 0.4, 'ticks': None, 'ticklabels': [] }
else:
	DEFAULT_KWARGS = {'m': wmap, 'plev': '800', 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
		'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 'hook': None,
		'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroscale': scale_oro,
		'oroalpha': 0.4, 'ticks': None, 'ticklabels': [] }
DEFAULT_CONTOUR_KWARGS = {'colors': 'k', 'alpha': 1.0, 'cmap': None, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'neither', 'linewidths': 2.0, 'linestyles': None }
DEFAULT_CONTOURF_KWARGS = {'colors': None, 'alpha': 1.0, 'cmap': None, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'neither', 'hatches': None }

MUTEX_GROUPS = [set(['colors', 'cmap']), ]

# DEFAULT settings per quantity Q
DEFAULT_Q = {}
DEFAULT_Q['defabs'] = {'cmap': _get_defabs_cm2(), 'extend': 'max', 'scale': scale_defabs, 'hook': hooks['defabs']}
DEFAULT_Q['defang'] = {'cmap': _get_periodic_cm3(), 'scale': scale_defang, 'ticks': ticks_defang, 
	'ticklabels': labls_defang}
DEFAULT_Q['u']   = {'scale': scale_u_diff}
DEFAULT_Q['Z']   = {'scale': scale_Z_diff}
DEFAULT_Q['T']   = {'cmap': plt.cm.RdBu_r,' extend': 'both'}
DEFAULT_Q['the']   = {'cmap': plt.cm.RdBu_r, 'extend': 'both'}
DEFAULT_Q['thestir']  = {'extend': 'both', 'cmap': plt.cm.BrBG}
DEFAULT_Q['thefil']  = {'extend': 'both', 'cmap': plt.cm.PRGn}
DEFAULT_Q['pv']  = {'scale': scale_pv, 'hook': hooks['pv']}
DEFAULT_Q['oro'] = {'scale': scale_oro_full, 'cmap': plt.cm.gist_earth, 'hook': hooks['oro']}
DEFAULT_Q['ow']  = {'scale': scale_ow_mean, 'extend': 'both'}
DEFAULT_Q['pvstir']  = {'scale': scale_pvstir_mean, 'extend': 'both'}
DEFAULT_Q['pvfil']  = {'scale': scale_pvfil_mean, 'extend': 'both', 'cmap': plt.cm.PRGn}
DEFAULT_Q['rsr']  = {'scale': scale_rsr_mean, 'extend': 'both'}
DEFAULT_Q['q'] = {'scale': scale_q, 'extend': 'max', 'cmap': _get_q_cm()}
DEFAULT_Q['qfil'] = {'scale': scale_qfs, 'extend': 'both', 'cmap': plt.cm.RdBu}
DEFAULT_Q['qstir'] = {'scale': scale_qfs, 'extend': 'both', 'cmap': plt.cm.BrBG}


# #############################################################################
# 5. Making the settings easily available
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
	default_q = DEFAULT_Q

	_overrides = {}

	def __init__(self):
		for q in Q:
			if q in self.default_q:
				self._overrides[q] = settings_dict(self.default_q[q], self)
			else:
				self._overrides[q] = settings_dict({}, self)

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
	

	def merge(self, q, **kwargs):
		rkwargs = dict(self.__getattribute__(q))
		for kwarg, argv in kwargs.items():
			rkwargs[kwarg] = argv

		return rkwargs


	def reset(self, q, key=None):
		self.__getattribute__(q).reset(key)

		return



class settings_contourf(settings_contour):
	default = settings_dict(DEFAULT_CONTOURF_KWARGS)
	default.update(DEFAULT_KWARGS)
	_overrides = {}
	


class settings(object):
	__default = {
		'q': Q,
		'bins': BINS_Q,
		'datapath': DATAPATH,
		'opath': OPATH,
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
		if key[0] == '_' or key in ['reset', ]:
			return object.__getattribute__(self, key)

		return self.__current[key]


	def __setattr__(self, key, value):
		if key in ['_settings__default', '__setattr__']:
			raise AttributeError, 'The default values cannot be overwritten'
		elif key[0] == '_' or key in ['reset', ]:
			return object.__setattr__(self, key, value)
		self.__current[key] = value

		return


	def reset(self, key=None):
		if key and key in self.__current:
			self.__current[key] = copy(self.__default[key])
		else:
			self.__current = deepcopy(self.__default)

		return

conf = settings()



# #############################################################################
# 6. Clean-Up: Making the default settings only available through settings objects
# 
del Q, BINS_Q, DATAPATH, OPATH, FILE_STD, FILE_STAT, FILE_MSTAT, STD_SLICE, YEARS, PLEVS, PTLEVS, PVLEVS
del DEFAULT_KWARGS, DEFAULT_CONTOUR_KWARGS, DEFAULT_CONTOURF_KWARGS, DEFAULT_Q, MUTEX_GROUPS


# that's it
