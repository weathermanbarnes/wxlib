#!/usr/bin/env python
# -*- encoding: utf-8


''' Defines additional data structures for plot settings storage and convenient access

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

 #. ``nd_default_dict`` objects: Derived from ``default_dict``, it adds an abitrary
    number of dimensions, accessible by ``dat[key1,key2,...,keyN]``. Valid keys for
    each dimension are prescribed, and default values must be provided for each key
    in the last dimesion. Default values are identical for all ``dat[::,keyN]``, 
    and hence independent from anything but the last dimension.
 #. ``plot_settings_dict``: A variant of the ``nd_default_dict`` used for storing
    plot configuration. The number of dimensions is fixed to three. They are, in order: 
    vertical level, variable and plot configuration key. 
'''

import numpy as np
from copy import copy 

from .. import proj
from .. import settings as base

import collections


class nd_default_dict(base.default_dict):
    ''' n-dimensional dictionary with default values for valid keys in the last dimension

    The default value is independent from any dimension but the last. All dimensions 
    before the last are subsequently called "first dimensions".

    The set of valid keys for the last dimension is limited, and the set of valid keys
    will typically not change once the object initiated. Valid keys cannot be removed, but
    in can be added with the ``add_default()`` method if necessary.

    Parameters
    ----------
    ndims : int
        The number of dimensions keys in this dictionary are required to have.
    defaults : dict
        Definition of valid keys for the last dimension, and their default values.
    mutexes : list of sets
        Sets of mutually exclusive configuration keys. As soon as one is not None, all others 
        are automaticall set to None.
    '''

    def __init__(self, ndims, defaults, mutexes):
        for default in defaults.values():
            if not isinstance(default, collections.Hashable):
                raise ValueError('Defaults for plotconf must be immutable, but got variable of type %s, value: %s' % (type(default), str(default)))
        self._ndims = ndims
        base.default_dict.__init__(self, defaults, mutexes)     # dictionary of valid keys in the last dimension and their default values
    
    def __check_key(self, key):
        ''' Check if a given query key is valid '''

        if len(key) < self._ndims - 1:
            raise KeyError('Key must contain at least %d dimensions' % self._ndims-1)
        elif len(key) > self._ndims:
            raise KeyError('Key can at most contain %d dimensions' % self._ndims)
        
        if len(key) == self._ndims and not key[-1] in self._defaults:
            raise KeyError(key)
        
        return


    def __default_walk(self, key, recur=None):
        ''' Return a list of keys to be checked for generic overrides'''

        if recur == None:
            recur = self._ndims-2

        if recur > 0:
            keys = []
            keys.extend(self.__default_walk(key, recur-1))
            if not key[recur] == None:
                nkey = key[:recur] + (None, ) + key[recur+1:]
                keys.extend(self.__default_walk(nkey, recur-1)) 
        else:
            if key[0] == None:
                keys = [key, ]
            else:
                nkey = (None, ) + key[1:]
                keys = [key, nkey]

        return keys

    
    def __getitem__(self, key):
        ''' Implement access via the ``dat[key1,key2,...,keyN]`` syntax '''

        self.__check_key(key)
        
        if len(key) == self._ndims:
            for key in self.__default_walk(key):
                if key in self._:
                    return self._[key]
            
            if key[-1] in self._defaults:
                return self._defaults[key[-1]]
            else:
                raise KeyError(key)
        else: 
            ret = {}
            for skey in self._defaults:
                fkey = key + (skey, )
                ret[skey] = self[fkey]

            return ret
    
    def __setitem__(self, key, value):
        ''' Implement assignment via the ``dat[key1,key2,...,keyN] = value`` syntax '''

        self.__check_key(key)

        # Replace slice(None) by None for multiple assignment
        for n, skey in zip(range(len(key)), copy(key)):
            if skey == slice(None):
                key = key[:n] + (None,) + key[n+1:]
        
        # Got one specific key to override
        if len(key) == self._ndims:
            # Respect that some combination configuration are not meaningful
            if type(value) == type(None):
                for mutex in self._mutexes.get(key[-1], []):
                    mkey = key[:-1] + (mutex,)
                    self._[mkey] = None
            self._[key] = value
        
        # Last dimension missing, expecting dict to override several keys at once
        else: 
            if not type(value) == dict:
                raise TypeError('Dict required for multiple assignment, got `%s` instead.' % type(value))

            for lkey, lvalue in value.items():
                fullkey = key + (lkey,)
                if not lkey in self._defaults:
                    raise KeyError(fullkey)
                # Respect that some combination configuration are not meaningful
                if type(lvalue) == type(None):
                    for mutex in self._mutexes.get(lkey, []):
                        mkey = key + (mutex,)
                        self._[mkey] = None
                self._[fullkey] = lvalue
                

    def __delitem__(self, key):
        ''' Implement reset via the ``del dat[key1,key2,...,keyN]`` syntax '''

        self.__check_key(key)

        for n, skey in zip(range(len(key)), copy(key)):
            if skey == slice(None):
                key = key[:n] + (None,) + key[n+1:]
        
        if len(key) == self._ndims:
            if key in self._:
                del self._[key]
        else:
            for lkey in self._defaults:
                fullkey = key + (lkey,)
                if fullkey in self._:
                    del self._[fullkey]

    def add_default(self, key, value):
        ''' Add a new key and its default value for the last dimension '''

        if key in self._defaults:
            raise KeyError('Default value for %s exists already' % str(key))

        self._defaults[key] = value


class plot_settings_dict(nd_default_dict):
    ''' A version of the nd_default_dict, where the number of dimensions is fixed to 3

    The dimensions are, in order: vertical level, variable and plot configuration key.

    Parameters
    ----------
    defaults : dict
        Definition of valid keys and their default values.
    mutexes : list of sets
        Sets of mutually exclusive configuration keys. As soon as one is not None, all others 
        are automaticall set to None.
    '''

    def __init__(self, defaults, mutexes):
        nd_default_dict.__init__(self, 3, defaults, mutexes)
    
    def merge(self, plev, q, **kwargs):
        ''' Merge given plot settings with the stored plot settings for (plev,q)

        Parameters
        ----------
        plev : str
            Vertical level name.
        q : str
            Variable name.

        Keyword arguments
        -----------------
        plot arguments : all
            For a list of valid arguments refer to :ref:`plot configuration`. Keyword 
            arguments are applied both to the inset map and the main cross section.

        Returns
        -------
        dict
            Merged plotconf parameters.
        '''

        merged = self[plev,q]
        merged.update(kwargs)
        
        return merged


# Make sure that the plot defaults are immutable, 
# such that they cannot be changed in place
PLOT_DEFAULTS = {
    'alpha': 1.0, 
    'cb_axes': None,
    'cb_disable': False, 
    'cb_expand_fig_fraction': 0.10,
    'cb_expand_pad_fraction': 0.2,
    'cb_label': None,
    'cb_orientation': 'horizontal',
    'cb_tickspacing': 'proportional',
    'cb_shrink': 0.8,
    'coastcolor': 'k', 
    'coastwidth': 1.0,
    'coast_alpha': None,
    'fig_size': 'auto',
    'fig_dpi': 150,
    'gridcolor': 'k',
    'grid_alpha': 0.3,
    'grid_dashes': (None, None),
    'grid_linestyle': '-',
    'grid_latmax': 90,
    'landcolor': None,
    'land_alpha': 1.0,
    'meridians': (-180,-120,-60, 0, 60, 120),
    'parallels': (-60,-30, 0, 30, 60),
    'hatches': (),
    'hook': None,
    'm': proj.world, 
    'maskcolor': '0.25',
    'mark': None, 
    'markersize': 25, # only for overlay_dots, TODO: homogeneize with mark keyword
    'name': '',
    'name_prefix': '',
    'oroalpha': 0.4, 
    'orocolor': 'k', 
    'oroscale': (9810., 19620., 29430., 39240., 49050., 58860., 68670., 78480.),
    'overlays': (), 
    'save': '', 
    '_return': None,
    'section_hor_resolution': 25000.0,
    'scale': 'auto', 
    'scale_exceed_percentiles': (0.01, 0.99),
    'scale_intervals': (1,2,3,5,10,),
    'scale_intervals_periodic': True,
    'scale_target_steps': 7,
    'scale_symmetric_zero': False,
    'show': True, 
    'ticks': None, 
    'ticklabels': (), 
    'title': 'auto', 
    'Zdata': None,
    'zorder': 1,

    'lon': None, 
    'lat': None, 
    'norm': None, 
    'vmin': None,
    'vmax': None, 
    'levels': None, 
    'origin': None, 
    'extent': None, 
}
PLOTL_DEFAULTS = copy(PLOT_DEFAULTS) 

PLOTL_DEFAULTS.update({
    'cmap': None, 
    'colors': 'k', 
    'contour_labels': False,
    'contour_labels_fontsize': 12,
    'contour_labels_inline': True,
    'contour_labels_inline_spacing': 2,
    'contour_labels_format': '%1.1f',
    'contour_labels_scale': None,
    'extend': 'neither', 
    'linestyles': None, 
    'linewidths': 2.0, 

    'vector_space_interval': 15,
    'vector_disable_interpolation': False,
    'vector_space_numbers_xy': (30,20),

    'barbcolor': None,
    'flagcolor': None,
    'length': None,
    'sizes': None,
    'fill_empty': None,
    'quiver_length': None,
    'quiver_width': None,
})
PLOT_DEFAULTS.update({
    'cmap': 'viridis',
    'colors': None, 
    'extend': 'both', 
    'tile': False,
    'tri': False,
    'edgecolors': 'none',
    
    'name_prefix': '',
    'plotformat': 'pdf',
    'plotpath': '.',
})

MUTEX_GROUPS = [set(['colors', 'cmap']), ]


# Create standard plot settings objects
pconfl = plot_settings_dict(PLOTL_DEFAULTS, MUTEX_GROUPS) 
pconf = plot_settings_dict(PLOT_DEFAULTS, MUTEX_GROUPS)


# C'est la fin.
