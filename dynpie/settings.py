#!/usr/bin/env python
# -*- encoding: utf-8

from copy import copy, deepcopy
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


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
# 2. Some often used map projections
#

# (a) World map
def wmap():
	return Basemap(projection='robin',lon_0=0,resolution='c')
# (b) Northern polar centered map
def npmap():
	return Basemap(projection='npstere',boundinglat=15,lon_0=-50,resolution='l')
# (c) Southern polar centered map
def spmap():
	return Basemap(projection='spstere',boundinglat=-15,lon_0=0,resolution='l')
# (d) North-Atlantic map
def NAmap():
	return Basemap(projection='lcc', lat_0=55, lat_ts=55, lon_0=-30, resolution='l', 
			width=9000000, height=6000000)
# (e) North-Pacific map
def NPmap():
	return Basemap(projection='lcc', lat_0=50, lat_ts=50, lon_0=-180, resolution='l', 
			width=9000000, height=6000000)
# (f) Australia map
def Ausmap():
	return Basemap(projection='lcc', lat_0=-50, lat_ts=-50, lon_0=120, resolution='l', 
			width=9000000, height=6000000)



# #############################################################################
# 3. Scales, ticks and labels
# 
scale_oro = range(10000,80001,10000)
scale_oro_full = range(-19000,51000,2000)

scale_Z_diff = np.arange(-5250,5251,500)

scale_u     = np.arange(20,71,10)
scale_u_diff = np.arange(-30,31,5)

scale_defabs = np.arange(5.0,30.1,5.0)
scale_defabs_mean = np.arange(2.0,12.1,2.0)

scale_defang = (np.arange(-18,19)-0.5)*np.pi/36.0
scale_defang_coarse = np.arange(-4,5)*np.pi/8.0 - np.pi/72.0
ticks_defang = np.arange(-4,5)*3.1415926535/8.0 
labls_defang = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']



# #############################################################################
# 4. Default settings
#
Q = {'defabs': 'defabs', 'defang': 'defang', 
	'm': 'mont', 'p': 'pres', 'u': 'u', 'v': 'v', 'Z': 'z', 'oro': 'oro'}
DATAPATH = ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY']
FILE_STD   = 'ei.ans.%d.%s.%s'
FILE_STAT  = 'ei.ans.%d.%s.%s.stat'
FILE_MSTAT = 'ei.ans.stat.%s.%s'
STD_SLICE  = (slice(None), slice(None), slice(None))

# DEFAULT contour settings
if os.getenv('DYNLIB_PLOT_PRINT'):
	DEFAULT_KWARGS = {'m': wmap, 'plev': 800, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
		'overlays': [], 'disable_cb': True, 'show': False, 'save': '', 'title': '',
		'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroscale': scale_oro,
		'oroalpha': 0.4, 'ticks': None, 'ticklabels': [] }
else:
	DEFAULT_KWARGS = {'m': wmap, 'plev': 800, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
		'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '',
		'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroscale': scale_oro,
		'oroalpha': 0.4, 'ticks': None, 'ticklabels': [] }
DEFAULT_CONTOUR_KWARGS = {'colors': 'k', 'alpha': 1.0, 'cmap': None, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'neither', 'linewidths': 2.0, 'linestyles': None }
DEFAULT_CONTOURF_KWARGS = {'colors': None, 'alpha': 1.0, 'cmap': None, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'neither', 'hatches': None }

# DEFAULT settings per quantity Q
DEFAULT_Q = {}
DEFAULT_Q['defabs'] = {'cmap': _get_defabs_cm2(), 'extend': 'max', 'scale': scale_defabs}
DEFAULT_Q['defang'] = {'cmap': _get_periodic_cm3(), 'scale': scale_defang, 'ticks': ticks_defang, 
	'ticklabels': labls_defang}
DEFAULT_Q['u'] = {'scale': scale_u_diff}
DEFAULT_Q['Z'] = {'scale': scale_Z_diff}
DEFAULT_Q['oro'] = {'scale': scale_oro_full, 'cmap': plt.cm.gist_earth}



# #############################################################################
# 5. Default settings
#
hooks = {}
hooks['defabs'] = lambda defabs: defabs*1e5
#hooks['oro'] = lambda oro: oro[



# #############################################################################
# 6. Making the settings easily available
#
class settings_contour(object):
	default = DEFAULT_CONTOUR_KWARGS
	default.update(DEFAULT_KWARGS)
	_overrides = {}

	def __init__(self):
		for q in Q:
			self._overrides[q] = copy(self.default)
			for setting in DEFAULT_Q.get(q, []):
				self._overrides[q][setting] = DEFAULT_Q[q][setting]

		return


	def __getattribute__(self, key):
		if key[0] == '_':
			return object.__getattribute__(self, key)
		elif key not in self._overrides:
			return object.__getattribute__(self, key)
		
		return self._overrides[key]


	def __setattr__(self, key, value):
		if key in self._overrides or key == 'default':
			raise AttributeError, 'The attributes cannot be overwritten'
		object.__setattr__(self, key, value)

		return
	

	def merge(self, key, **kwargs):
		rkwargs = copy(self.__getattribute__(key))
		for kwarg, argv in kwargs.items():
			rkwargs[kwarg] = argv

		return rkwargs


	def reset(self, key, setting=None):
		if not setting:
			self._overrides[key] = copy(self.default)
		else:
			self._overrides[key][setting] = self.default[setting]

		return



class settings_contourf(settings_contour):
	default = DEFAULT_CONTOURF_KWARGS
	default.update(DEFAULT_KWARGS)
	_overrides = {}
	


class settings(object):
	__default = {
		'q': Q,
		'datapath': DATAPATH,
		'file_std': FILE_STD,
		'file_stat': FILE_STAT,
		'file_mstat': FILE_MSTAT,
		'std_slice': STD_SLICE,
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
# 7. Clean-Up: Making the default settings only available through settings objects
# 
del Q, DATAPATH, FILE_STD, FILE_STAT, FILE_MSTAT, STD_SLICE
del DEFAULT_KWARGS, DEFAULT_CONTOUR_KWARGS, DEFAULT_CONTOURF_KWARGS, DEFAULT_Q


# that's it
