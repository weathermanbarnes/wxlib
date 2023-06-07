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

import numpy as np
from copy import copy 
from .metio.standard_variables import standard_variables

import collections
import sys
if sys.version_info.major == 3 and sys.version_info.minor >= 10:
    from collections.abc import MutableMapping, Hashable
else:
    from collections import MutableMapping, Hashable


class default_dict(MutableMapping):
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
            if not isinstance(value, Hashable):
                self._[key] = copy(value)

        MutableMapping.__init__(self)
        
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
            if not isinstance(self._defaults[key], Hashable):
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
    disabled and provides identical functionality.

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
    _Q_LINES = 'q_lines'
    _Q_OBJ = 'q_obj'
    _Q_STD = 'q_std'
    _Q_AVG = 'q_avg'
    
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
    
    def __register_single_variable(self, var):
        ''' Add a single variable to the configurarion '''
        
        q = var.q

        # Long name and units must be set
        self[self._Q_LONG][q] = var.long_name
        self[self._Q_UNITS][q] = var.units

        # Optional: file segment
        if not type(var.qf) == type(None) and not var.qf == q:
            self[self._Q_FILE][q] = var.qf
            self[self._Q][var.qf] = q

        # Optional: mapping to a standard shorthand 
        if not type(var.q_std) == type(None) and not var.q_std == q:
            self[self._Q_STD][q] = var.q_std

        # Optional: mapping to a time-averge name
        if not type(var.q_avg) == type(None) and not var.q_avg == q:
            self[self._Q_AVG][q] = var.q_avg

        # Optional: Register as special variable
        if not type(var.bins) == type(None):
            self[self._Q_BINS][q] = var.bins
        if not type(var.lines) == type(None):
            self[self._Q_LINES][q] = var.lines
            self[self._Q_FILE][var.lines] = q
        if not type(var.objmask) == type(None):
            self[self._Q_OBJ][q] = var.objmask


    def register_variable(self, variables):
        ''' Register (a) new variable(s) 
        
        Parameters
        ----------
        variables : (list of) variable object(s)
            A single or list of variable definitions, using the variable object defined in the generic data source.
        '''

        if type(variables) == list:
            for var in variables:
                self.__register_single_variable(var)
        else:
            self.__register_single_variable(variables)
        

# Not every key is applicable and used for all data sources -- but no additional keys than 
# those listed here should be defined per data source.
default_conf = settings_obj({
    'q': {},        # Given a file name segment, which variable to we expect to find there?
    'qf': {},       # Given a variable name, in which file do we find it?
    'q_std': {},     # Given a data source-specific variable name, which standard variable name does it correspond to? (used primarily for plotting currently)
    'q_avg': {},    # Given a variable name, what is the name of a time-averge (if different)
    'q_units': {},
    'q_long': {},
    'q_bins': {},
    'q_lines': {},
    'q_feature_dists': {},
    'q_obj': {},
    'datapath': ['.', '/Data/gfi/users/local/share'], 
    'staticfile': '',
    'opath': '.',
    'oformat': 'nc',
    'epoch': None,
    'timestep': None,
    'gridsize': None,
    'local_timezone': 'Europe/Oslo',
}, [])
default_conf.register_variable(standard_variables)


# C'est la fin.
