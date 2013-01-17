#!/usr/bin/env python
# -*- encoding: utf-8

from copy import copy, deepcopy


# Default settings
Q = {'defabs': 'defabs', 'defang': 'defang', 
	'm': 'mont', 'p': 'pres', 'u': 'u', 'v': 'v', 'Z': 'z' }
DATAPATH = ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY']
FILE_STD   = 'ei.ans.%d.%s.%s'
FILE_STAT  = 'ei.ans.%d.%s.%s.stat'
FILE_MSTAT = 'ei.ans.stat.%s.%s'
STD_SLICE  = slice(None)

# DEFAULT contour settings
DEFAULT_KWARGS = {'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
	'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '',
	'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k',
	'oroalpha': 0.4, 'ticks': None, 'ticklabels': [] }
DEFAULT_CONTOUR_KWARGS = {'colors': None, 'alpha': 1.0, 'cmap': None, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'neither', 'linewidths': None, 'linestyles': None }
DEFAULT_CONTOURF_KWARGS = {'colors': None, 'alpha': 1.0, 'cmap': None, 'norm': None, 
	'vmin': None, 'vmax': None, 'levels': None, 'origin': None, 'extent': None, 
	'extend': 'neither', 'hatches': None }



class settings_contour(object):
	default = DEFAULT_CONTOUR_KWARGS
	default.update(DEFAULT_KWARGS)
	_overrides = {}

	def __init__(self):
		for q in Q:
			self._overrides[q] = copy(self.default)

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

s = settings()



# Making the default settings only available through settings objects
del Q, DATAPATH, FILE_STD, FILE_STAT, FILE_MSTAT, STD_SLICE
del DEFAULT_KWARGS


# that's it
