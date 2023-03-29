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

from ... import utils



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

    def __init__(self, name, plev=None, q=None):
        self.name = name
        self.requires = None
        
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
        if type(self.requires) == type(None):
            if type(b.requires) == type(None):
                ret.get_time_series = lambda dates: np.logical_and(self.get_time_series(dates), b.get_time_series(dates))
                ret.requires = None
            else:
                ret.get_time_series = lambda dates, *stuff: np.logical_and(self.get_time_series(dates), 
                        b.get_time_series(dates, *stuff))
                ret.requires = b.requires
        else:
            if type(b.requires) == type(None):
                ret.get_time_series = lambda dates, *stuff: np.logical_and(self.get_time_series(dates, *stuff), 
                        b.get_time_series(dates))
                ret.requires = self.requires
            else:
                ret.get_time_series = lambda dates, *stuff: np.logical_and(self.get_time_series(dates, *stuff), 
                        b.get_time_series(dates, *stuff))
                ret.requires = self.requires + b.requires

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
            True if the respective date should be part of the composite, and False otherwise.
        '''

        raise NotImplementedError('`decider.get_time_series` must be overriden in derived classes!')



class decide_by_data(decider):
    ''' Compositing criteria based on a given time series

    Parameters
    ----------
    name : str
        Name of the composite, to be used for saving.
    q : str
        The variable name identifier for required test data, following the data source
        conventions, e.g. ``'u'`` or ``'msl'`` in ERA-Interim.
    plev : str
        The vertical level of the required test data, following the ECMWF conventions, i
        e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface.
    criterion : function
        A function mapping a 3d-snapshot of q at plev(s) to True/False.
    '''

    def __init__(self, name, plev, q, criterion):
        super().__init__(name, plev, q)
        self.requires = (plev, q)
        self.evaluate = criterion

        return

    def get_time_series(self, dates, files_by_plevq, get_from_file, kwargs):
        __doc__ = super().get_time_series.__doc__
        
        start, end = min(dates), max(dates)+td(0,1) # Final date should be included in the request below
        req = list(files_by_plevq((self.plev, self.q), start=start, end=end))
        datshape = req[0][3][1:]      # Shape of the resulting data arrays
        for entry in req[1:]:
            shape = entry[3]
            if not shape[1:] == datshape[plevq]:
                raise ValueError(f'''Discovered inconsistent data shape across time:
                        plevq: {plevq}
                        file {entry[0]} with shape {shape[1:]}, 
                        preceeding files with shape {datshape[plevq]}.''')

        comp_ts = []
        for filename, tidxs, dates_, shape in req:
            cut = slice(tidxs[0], tidxs[-1]+1)
            tlen = len(tidxs)
            
            dat_ = get_from_file(filename, self.plev, self.q, cut=cut, no_static=True, **kwargs)
            
            # TODO: Treat lines 
            #if self.q in LINES:
            #    datoff_, grid = get_from_file(filename, self.plev, LINES[self.q], cut=cut)
            #    dat_ = utils.normalize_lines(dat_, datoff_, grid.dx, grid.dy)[:,np.newaxis,:,:]

            # Object ID masks are kept as is
            # Binned variables are kept as they are

            for tidx, date in enumerate(dates_):
                if date not in dates:
                    continue

                comp_ts.append(self.evaluate(dat_[tidx,::]))

        if not len(comp_ts) == len(dates):
            raise ValueError(f'{len(dates)-len(comp_ts)} time steps of the requested dates is not available for the test data.')

        return np.array(comp_ts)



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
        self.tidx = 0

        return

    def get_time_series(self, dates):
        __doc__ = super().get_time_series.__doc__
        
        comp_ts = []
        for date in dates:
            # Out of bounds -> always return False
            if date < self.dates[0] or date >= self.dates[-1]:
                comp_ts.append(False)

            else:
                # Match date with the largest one in self.dates that is smaller than date
                while self.dates[self.tidx+1] <= date:
                    self.tidx += 1

                # Evaluate the corresponding value
                comp_ts.append(self.decisions[self.tidx])

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
    tl_decider.name += f'_lag{sgn}{int(abs(tdiff.total_seconds())/3600):03d}h'

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

    tl_deciders = [__timelag_one(decider, tdiff) for tdiff in tdiffs]

    return tl_deciders



def __timelag_one_relativedelta(decider, factor, unit_interval):
    ''' Helper function to create a variant of a decider with time lag specified through relativedeltas '''
    
    tdiff = factor * unit_interval
    tl_decider = deepcopy(decider)
    tl_decider.__get_time_series = tl_decider.get_time_series
    
    def new_get_time_series(dates):
        # cftime objects need to be converted to datetime for dateutils to work
        if hasattr(dates[0], '_to_real_datetime'):
            dates = [d._to_real_datetime() for d in dates]
        return tl_decider.__get_time_series([date-tdiff for date in dates])

    tl_decider.get_time_series = new_get_time_series

    if factor >= 0:
        sgn = '+'
    else: 
        sgn = '-'

    if unit_interval.months > 0:
        number = abs(factor)*unit_interval.months
        unit = 'm'
    elif unit_interval.years > 0:
        number = abs(factor)*unit_interval.years
        unit = 'y'
    else:
        raise ValueError('Unit must be a positive number of years or months.')

    tl_decider.name += f'_lag{sgn}{number:03d}{unit}'

    return tl_decider



def timelag_relativedelta(decider, factors, unit_interval):
    ''' Create a series of time-lagged versions of one decider based on relativedelta objects

    Parameters
    ----------
    decider : decider
        Decider to be time-lagged.
    factors : list of int
        Factors by which the unit_interval is to be multiplied to arrive at the time lag.
    unit_interval : dateutil.relativedelta
        Unit time interval (e.g. 1 month) defining the lags.

    Returns
    -------
    dict of list of decider
        List of time-lagged versions of the decider, wrapped in a named group
        such that they will by default be saved together.
    '''

    tl_deciders = [__timelag_one_relativedelta(decider, factor, unit_interval) for factor in factors]

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
    
    return [item1 & item2 for item1 in list1 for item2 in list2]



# C'est le fin
