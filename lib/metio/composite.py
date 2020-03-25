#!/usr/bin/env python
# -*- encoding: utf-8


''' This module defines decider functions for composites

Each decider implements a criterion to deterime if a given time step should 
be part of the composite or not. The decision can be based on either time
itself, a provided time series or test data.

Several deciders can be combined by logical operators. For example a combination
of a decider for positive NAO phases and another decider for winter can be used
to construct a combined decider for positive NAO during winter.

Furthermore, the module provides functions to create lagged composites from a
given decider, and to create all combinations between two lists of composites
(e.g. for combining all variability indexes with all seasons).
'''


import numpy as np

from copy import deepcopy
from datetime import timedelta as td

from ..settings_basic import conf



# General base class
class decider(object):
    ''' Abstract base class for all deciders
    
    The base class defines the shared API for deciders implements the logical 
    operators to combine deciders. 

    Parameters
    ----------
    name : str
        Name of the composite, to be used for saving.
    q : str
        *Optional*. If applicable, the variable name identifier for required test 
        data, following the ECMWF conventions, e.g. ``'u'`` or ``'msl'``. 
    plev : str
        *Optional*. If applicable, the vertical level of the required test data, 
        following the ECMWF conventions, e.g. ``'700'`` for 700 hPa or ``'pv2000'`` 
        for the PV2-surface.
    '''

    def __init__(self, name, q=None, plev=None):
        self.name = name
        
        if q and plev:
            self.required_qs = (plev, q)

        if plev:
            self.plev = plev
            self.name += '_%s' % plev
        if q:
            self.q = q
            self.name += '_%s' % q

        return

    
    def __invert__(self):
        ''' Implements the ``~deciderA`` syntax '''

        ret = decider(f'!{self.name}')
        ret.get_time_series = lambda dates: np.logical_not(self.get_time_series(dates))

        return ret

    def __xor__(self, b):
        ''' Implements the ``deciderA ^ deciderB`` syntax '''

        ret = decider(f'{self.name}_xor_{b.name}')
        ret.get_time_series = lambda dates: np.logical_xor(self.get_time_series(dates), b.get_time_series(dates))

        return ret

    def __or__(self, b):
        ''' Implements the ``deciderA | deciderB`` syntax '''

        ret = decider(self.name+b.name)
        ret.get_time_series = lambda dates: np.logical_or(self.get_time_series(dates), b.get_time_series(dates))

        return ret

    def __and__(self, b):
        ''' Implements the ``deciderA & deciderB`` syntax '''

        ret = decider(self.name+'@'+b.name)
        ret.get_time_series = lambda dates: np.logical_and(self.get_time_series(dates), b.get_time_series(dates))

        return ret
    
    # The deciding function, returns a boolean time series for the dates given
    def get_time_series(self, dates):
        ''' The actual deciding function if a given time should be part of the composite

        Parameters
        ----------
        dates : list of datetime
            Date of the given time
    
        Returns 
        -------
        np.ndarray of dtype bool
            True of the respective date should be part of the composite, and False otherwise.
        '''

        raise NotImplementedError('`decider.get_time_series` must be overriden in derived classes!')



class decide_by_timeseries(decider):
    ''' Compositing criteria based on a given time series

    Parameters
    ----------
    name : str
        Name of the composite, to be used for saving.
    ts : dict
        Time series dictionary, containing a list of datetime objects and a list of values.
    criterion : function
        A function mapping a time series value to True/False.
    '''

    def __init__(self, name, ts, criterion):
        super().__init__(name)
        self.dates  = ts['dates']
        self.decisions = criterion(ts['values'])

        return

    def get_time_series(self, dates):
        __doc__ = super().get_time_series.__doc__
        
        comp_ts = []
        tidx = 0
        for date in dates:
            # Out of bounds -> always return False
            if date < self.dates[0] or date >= self.dates[-1]:
                comp_ts.append(False)

            else:
                # Match date with the largest one in self.dates that is smaller than date
                while self.dates[tidx+1] <= date:
                    tidx += 1

                # Evaluate the corresponding value
                comp_ts.append(self.decisions[tidx])

        return np.array(comp_ts)



class decide_by_date(decider):
    ''' Compositing criteria based on the date ifself

    Parameters
    ----------
    name : str
        Name of the composite, to be used for saving.
    criterion : function
        A function mapping a date to True/False.
    '''

    def __init__(self, name, criterion):
        super().__init__(name)
        self.evaluate = criterion

        return

    def get_time_series(self, dates):
        __doc__ = super().get_time_series.__doc__
        
        return np.array([self.evaluate(date) for date in dates])



def __timelag_one(decider, tdiff):
    ''' Helper function to create a variant of a decider with lag tdiff '''

    tl_decider = deepcopy(decider)
    tl_decider.__get_time_series = tl_decider.get_time_series
    tl_decider.get_time_series = lambda dates: tl_decider.__get_time_series(
            [date-tdiff for date in dates])

    if tdiff >= td(0): 
        sgn = '+'
    else: 
        sgn = '-'
    tl_decider.name += f'{sgn}{int(abs(tdiff.total_seconds())/3600):03d}h'

    return tl_decider



def timelag(decider, tdiffs):
    ''' Create a series of time-lagged versions of one decider 

    Parameters
    ----------
    decider : decider
        Decider to be time-lagged.
    tdiffs : list of timedelta
        List of timedelta-objects by which the decider should be time-lagged.

    Returns
    -------
    dict of list of decider
        List of time-lagged versions of the decider, wrapped in a named group
        such that they will by default be saved together.
    '''

    tl_deciders = {decider.name: [__timelag_one(decider, tdiff) for tdiff in tdiffs]}

    return tl_deciders



def matrix(list1, list2):
    ''' Create a list of deciders containing all combinations of the given decider lists
    
    If the two lists are, for example, a list of climate variability indexes and a list 
    of seasons, a list of deciders for all variability indexes for all seasons is returned.

    Parameters
    ----------
    list1: list of decider
        First list of deciders to be combined.
    list2: list of decider
        Second list of deciders to be combined. Each decider in this list makes a resulting category.

    Returns
    -------
    list of decider
        List of combinations.
    '''

    # Flatten dicts
    if type(list1) == dict:
        list1 = [item for group in list1.values() for item in group]
    if type(list2) == dict:
        list2 = [item for group in list2.values() for item in group]
    
    # Wrap single deciders in a list
    if not type(list1) == list: list1 = [list1, ]
    if not type(list2) == list: list2 = [list2, ]
    

    if len(list1) == 0 or len(list2) == 0:
        return {}
    
    return {item2.name: [item1 & item2 for item1 in list1] for item2 in list2}



# C'est le fin
