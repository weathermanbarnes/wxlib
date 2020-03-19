#!/usr/bin/env python
# -*- encoding: utf-8


import os
import netCDF4 as nc
import scipy.io.matlab as mat

import numpy as np
import cftime
from datetime import timedelta as td

from ..settings_basic import conf
from ..gridlib import grid_by_static, grid_by_nc
from .. import utils


# 16 billion values -> 128G of memory, a bit more than 40 years of 6 hourly data at 0.5deg resolution
MAX_REQUEST_SIZE = 16.0e9 
WARN_REQUEST_SIZE = 2.0e9


timestep = None
gridsize = None
staticfile = None
calendar = 'standard'
dt = cftime.DatetimeGregorian



def plev_to_tuple(plev, with_unit=False):
    ''' Return meta-information about a given plev string

    Returns the type of vertical level and a numeric value

    Parameters
    ----------
    plev : str
        The level specifying string.
    with_unit: bool
        If ``'True'``, return string giving the unit as a third return variable.

    Returns
    -------
    plevtype : str
        The type of vertical level, e.g. ``'isobaric'``, ``'isentropic'``, 
        ``'potential vorticity'``, ``'height above ground'``.
    plevhgt : float
        The height of the vertical level in the appropriate unit for the ``'plevtype'``.
    plevunit : str
        If ``'with_unit'``, return unit for the height of the vertical level.
    '''

    if plev[:2] == 'pt':
        plevtype = 'isentropic'
        plevunit = 'K'
        plevhgt = float(plev[2:])
    elif plev[:2] == 'pv':
        plevtype = 'potential vorticity'
        plevunit = 'PVU'
        plevhgt = float(plev[2:])/1000.0
    elif plev[-1:] == 'm':
        plevtype = 'height above ground'
        plevunit = 'm'
        plevhgt = float(plev[:-1])
    else:
        try: 
            plevhgt = float(plev) * 100.0
        except ValueError:
            raise Exception(f'Cannot interpret vertical level string `{plev}`')
        plevtype = 'isobaric'
        plevunit = 'Pa'

    if with_unit:
        return plevtype, plevhgt, plevunit
    else:
        return plevtype, plevhgt


def tuple_to_plev(plevtype, plevhgt):
    ''' Return brief plev descriptor string based on type and height

    Parameters
    ----------
    plevtype : str
        The type of vertical level, e.g. ``'isobaric'``, ``'isentropic'``, 
        ``'potential vorticity'``, ``'height above ground'``.
    plevhgt : float
        The height of the vertical level in the appropriate unit for the ``'plevtype'``.

    Returns
    -------
    plev : str
        The level specifying string.
    '''

    if plevtype == 'isobaric':
        plev = str(int(plevhgt/100))
    elif plevtype == 'isentropic':
        plev = f'pt{int(plevhgt)}'
    elif plevtype == 'potential vorticity':
        plev = f'pv{int(plevhgt*1000)}'
    elif plevtype == 'height above ground':
        plev = f'{int(plevhgt)}m'
    else:
        raise NotImplementedError(f'Unknown vertical level type `{plevtype}`')

    return plev


class files_by_plevq(object):
    def __init__(self, plevq, start, end):
        plev, q = plevq
        self.plev = plev
        self.q = q
        self.start = start
        self.cur = start
        self.end = end

        return 

    def __iter__(self, ):
        return self
    
    def __next__(self):
        raise StopIteration


# TODO: Return xarray by default, reflect other docu changes (if any?)
def metopen_factory(get_static):
    ''' Create the metopen function based on data source specific helpers '''

    def handle_npy(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read npy files '''

        if not mode == 'r':
            print('WARNING: Can only open npy files in read mode!')
        if q:
            dat = np.load(filepath, mmap_mode='r', encoding='bytes')
            dat = dat[cut]
        else:
            dat = None
        if not quiet:
            print('Found '+filepath)
        f = None

        return f, dat

    def handle_npz(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read npz files '''

        if not mode == 'r' and not quiet:
            print('WARNING: Can only open npz files in read mode!')
        f = np.load(filepath, encoding='bytes')
        if q:
            if q not in f.files:
                tried.append(filepath)
            dat = f[q][cut]
        else:
            dat = None
        if not quiet:
            print('Found '+filepath)

        return f, dat
    
    def handle_mat(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read mat files '''

        if not mode == 'r' and not quiet:
            print('WARNING: Can only open mat files in read mode!')
        f = mat.loadmat(filepath)
        if q:
            if q not in f:
                tried.append(filepath)
            dat = f[q][cut]
        else:
            dat = None
        if not quiet:
            print('Found '+filepath)

        return f, dat

    def handle_nc(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read netCDF files '''

        f = nc.Dataset(filepath, mode)
        f.set_auto_scale(False)
        f.set_auto_mask(False)
        if q:
            var = f.variables[q]
            if q not in f.variables:
                tried.append(filepath)
            if not no_dtype_conversion:
                dat = utils.scale(var, cut=cut)
            else:
                dat = var[cut]
        else:
            dat = None

        if not no_static:
            if q:
                static = grid_by_nc(f, var)
            else:
                static = grid_by_nc(f)
            # TODO: Where to search for topography in nc files?
            static.oro = np.zeros((static.ny, static.nx))
        else:
            static = None

        if not quiet:
            print('Found '+filepath)

        return f, dat, static


    def metopen(filename, q=None, **kwargs):
        ''' Find and open files by name
        
        Uses the conf.datapath list to locate files in a variety of different locations. 
        In user scripts and user settings files, this variable will typically be available
        via conf.datapath.

        The located files might either be numpy-files, netCDF-files or matlab mat-files. 
        For each of the file types, metopen returns the requested variable and some meta-
        information about the variable, if not suppressed.

        Parameters
        ----------
        filename : str
            The name of the file, excluding the file ending.
        q : str
            *Optional*. The requested variable within the file.
        cut : slice
            *Optional*, default ``slice(None)``. Limit the request to a given time slice. 
            With the default data layout, only relevant data needs to be read when only 
            a time slice of the entire data is requested. Hence, using cut to limit your 
            data request can make reading the data largely more efficient.
        verbose : bool
            *Optional*, default ``False``. Print debug information on which files are 
            being looked for.
        no_dtype_conversion : bool
            *Optional*, default ``False``. By default, ``metopen`` uncompresses data in the 
            file and converts all data to float64. This behaviour can be suppressed by 
            setting this option to ``True``. This implies, however, that scaling and offset
            cannot be applied automatically.
        no_static : bool
            *Optional*, default ``False``. By default, ``metopen`` does its best to 
            provide meta-information about the requested file, using 
            :module:`grid.gridlib`, and returns the meta-information as a third value. 
            This behaviour can be suppressed by setting this parameter to ``True``.
        mode : str
            *Optional*, default ``'r'``. Only effective for netCDF files. The read/write mode
            with which to open the file. Valid values are ``'r'`` for read-only access, `'a'``
            or ``'r+'` for read-write access and and ``'w'``for replacing the given file. 
        no_xarray : bool
            *Optional*, default ``False``. If mode='r', netCDF files are by default returned 
            as xarray data sets. If ``True`` metopen instead returns an netCDF4 file object.

        Returns
        -------
        data file object
            xarray dataset or python data, netCDF or Matlab file object.
        np.ndarray
            If q given, data of the requested variable.
        grid.gridlib
            If ``no_static=False`` meta-information about the requested data.
        '''

        kwargs['no_static'] = kwargs.get('no_static', False)
        kwargs['verbose'] = kwargs.get('verbose', False)
        kwargs['quiet'] = kwargs.get('quiet', False)
        kwargs['no_dtype_conversion'] = kwargs.get('no_dtype_conversion', False)
        kwargs['cut'] = kwargs.get('cut', slice(None))
        
        if not type(kwargs['cut']) == slice:
            raise ValueError('cut must be a 1-dimensional slice object')

        tried = []
        for path in conf.datapath:
            static = None

            if kwargs['verbose']:
                print('Trying: '+path+'/'+filename+'.*')

            if os.path.exists(path+'/'+filename+'.npy'):
                f, dat = handle_npy(path+'/'+filename+'.npy', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'npy':
                f, dat = handle_npy(path+'/'+filename, q, **kwargs)

            elif os.path.exists(path+'/'+filename+'.npz'):
                f, dat = handle_npz(path+'/'+filename+'.npz', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'npz':
                f, dat = handle_npz(path+'/'+filename, q, **kwargs)

            elif os.path.exists(path+'/'+filename+'.mat'):
                f, dat = handle_mat(path+'/'+filename+'.mat', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'mat':
                f, dat = handle_mat(path+'/'+filename, q, **kwargs)

            elif os.path.exists(path+'/'+filename+'.nc'):
                f, dat, static = handle_nc(path+'/'+filename+'.nc', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'nc':
                f, dat, static = handle_nc(path+'/'+filename, q, **kwargs)

            else:
                tried.append(path)
                continue
            
            if q and not kwargs['no_dtype_conversion'] and not dat.dtype == 'f8':
                dat = dat.astype('f8')

            if not kwargs['no_static']:
                if not static:
                    static = get_static(kwargs['verbose'], kwargs['no_dtype_conversion'], kwargs['quiet'])
            else:
                if q:
                    return f, dat
                else:
                    return f
            if q:
                return f, dat, static
            else:
                return f, static
        
        raise ValueError('%s.* not found in any data location. \n'
                'Tried the following (in order):\n\t%s' % (filename, '\n\t'.join(tried)) )

    return metopen



def metsave(**kwargs):
    # TODO: Should may be become a generic saving function similar to xarray.to_netcdf(), 
    # but with proper group support?
    #
    # In every case: Warn about missing meta data rather than stopping with errors.
    ''' TBD '''

    raise NotImplementedError('TBD. May be data source-specific, may be not.')



def get_instantaneous_factory(metopen, files_by_plevq, get_static):
    ''' Create the get_instantaneous function based on data source specific helpers '''

    def get_instantaneous(plevqs, dates, force=False, **kwargs):
        ''' Get instantaneous fields

        Allows general data requests in the configured data base, e.g. ERA-Interim. The request
        can span several files, e.g. by including several vertical levels or by covering several
        years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
        
        The method internally uses metopen to locate data files. Hence, it will find data in the 
        locations given conf.datapath (in user scripts and user settings files typically 
        available as conf.datapath).

        Parameters
        ----------
        plevqs : 2-tuple or list of 2-tuples
            Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
            e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
            identifier, following the ECMWF conventions as far as appicable, e.g. ``'u'`` or ``'msl'``.
            Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
            and/or the variable name, to request all vertical levels/variables available.
        dates : list of datetime
            The minimum and maxmimum dates in this list define the requested time interval. The i
            final date will not be included in the result, i.e. for all time steps in 2016 request
            dates from 2016-01-01 00:00 to 2017-01-01 00:00.
        force : bool
            *Optional*, default ``False``. Turn off the error, if large amounts of data are
            requested at once. **Be sure you know what you are doing, when setting this to 
            ``True``! Your request might make your script occupy a large fraction of the 
            system memory**.
        
        Keyword arguments
        -----------------
        metopen arguments : all optional
            Optional arguments passed on to calls of metopen within this function.

        Returns
        -------
        dict of np.ndarray
            Data for the requested variable.
        grid.gridlib
            If ``no_static=False`` meta-information about the requested data, otherwise ``None``.
        '''
        
        start, end = min(dates), max(dates)

        if type(plevqs) == tuple:
            plevqs = [plevqs, ]

        # Dry run to test whether input parameters are valid and to determine the total amount of requested data
        req = {}            # Request info
        datshape = {}       # Shape of the resulting data arrays
        request_size = 0    # Total size of all requested data 
        for plevq in plevqs:
            # req[plevq] contains a list of 4-tuples: 
            # (filename, list of tidx, list of dates, request_size in number of values)
            req[plevq] = list(files_by_plevq(plevq, start=start, end=end))
            datshape[plevq] = req[plevq][0][3]
            for entry in req[plevq][1:]:
                shape = entry[3]
                if not shape[1:] == datshape[plevq][1:]:
                    raise ValueError(f'''Discovered inconsistent data shape across time:
                            plevq: {plevq}
                            file {entry[0]} with shape {shape[1:]}, 
                            preceeding files with shape {datshape[plevq][1:]}.''')
                datshape[plevq] = (datshape[plevq][0] + shape[0],) + shape[1:]
            request_size += sum([np.prod(reqinfo[3]) for reqinfo in req[plevq]])

        # Checking max length
        if request_size > MAX_REQUEST_SIZE:
            if force:
                print(f'Warning: you requested {request_size / 1024**3:.2f}G of data.')
            else:
                raise ValueError('''Cowardly refusing to fetch {request_size} values at once.
                    If you are absolutely certain what you're doing, you can use force=True to override.''')

        if not force and request_size > WARN_REQUEST_SIZE:
            print(f'Warning: you requested {request_size / 1024**3:.2f}G of data.')
        
        # Remove no_static if present, this should only effect the output of this function rather than the metopen calls herein
        no_static = kwargs.pop('no_static', False)

        dat = {}
        dates = {}
        for plev, q in plevqs:
            # Allocate memory for the result
            dat[plev,q] = np.empty(datshape[plev,q])
            dates[plev,q] = []
            
            toff = 0
            for filename, tidxs, dates_, shape in req[plev,q]:
                cut = slice(tidxs[0], tidxs[-1]+1)
                tlen = len(tidxs)
                
                # TODO: ERA5 will in general need an interpolation at this stage. 
                #       How to incorporate in the general mechanism?
                f, dat_ = metopen(filename, q, cut=cut, no_static=True, **kwargs)
                dat[plev,q][toff:toff+tlen,...] = dat_[...]
                dates[plev,q].extend(dates_)
                
                toff += tlen

        # TODO: Check whether dates lists are consistent (?)
        
        # Prepare grid information, if so requested
        if no_static:
            grid = None
        else:
            grid = get_static()
            grid = grid.new_time(dates[plev,q])
            
        return dat, grid
    
    return get_instantaneous


def get_time_average(plevqs, dates, **kwargs):
    ''' Get time-average fields

    Allows general data requests in the configured data base, e.g. ERA-Interim. The request
    can span several files, e.g. by including several vertical levels or by covering several
    years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
    
    The method internally uses metopen to locate data files. Hence, it will find data in the 
    locations given conf.datapath (in user scripts and user settings files typically 
    available as conf.datapath).

    Parameters
    ----------
    plevqs : 2-tuple or list of 2-tuples
        Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
        e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
        identifier, following the ECMWF conventions as far as appicable, e.g. ``'u'`` or ``'msl'``.
        Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
        and/or the variable name, to request all vertical levels/variables available.
    dates : list of datetime
        The minimum and maxmimum dates in this list define the requested time interval. The i
        final date will not be included in the result, i.e. for an average of 2016 request
        dates from 2016-01-01 00:00 to 2017-01-01 00:00.
    
    Keyword arguments
    -----------------
    metopen arguments : all optional
        Optional arguments passed on to calls of metopen within this function.

    Returns
    -------
    dict of np.ndarray
        Time average for the requested variable and vertical level(s).
    grid.gridlib
        If ``no_static=False`` meta-information about the requested data, otherwise ``None``.
    '''

    raise NotImplementedError('to be written')


def get_aggregate(plevqs, dates, agg, **kwargs):
    ''' Get time-aggregates, such as monthly means

    Allows general data requests in the configured data base, e.g. ERA-Interim. The request
    can span several files, e.g. by including several vertical levels or by covering several
    years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
    
    The method internally uses metopen to locate data files. Hence, it will find data in the 
    locations given conf.datapath (in user scripts and user settings files typically 
    available as conf.datapath).

    Parameters
    ----------
    plevqs : 2-tuple or list of 2-tuples
        Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
        e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
        identifier, following the ECMWF conventions as far as appicable, e.g. ``'u'`` or ``'msl'``.
        Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
        and/or the variable name, to request all vertical levels/variables available.
    dates : list of datetime
        The minimum and maxmimum dates in this list define the requested time interval. The i
        final date will not be included in the result, i.e. for all time steps in 2016 request
        dates from 2016-01-01 00:00 to 2017-01-01 00:00.
    agg : str or dynlib.tagg.agg object
        String representation (e.g. ``'cal_month'``, ``'pentad'``, or ``'met_season'``), or 
        time aggregator object representation the aggregation interval.
    
    Keyword arguments
    -----------------
    metopen arguments : all optional
        Optional arguments passed on to calls of metopen within this function.

    Returns
    -------
    dict of np.ndarray
        Time-aggregated data for the requested variable(s) and vertical level(s).
    grid.gridlib
        If ``no_static=False`` meta-information about the requested data, otherwise ``None``.
    '''

    raise NotImplementedError('to be written')


def get_composite(plevqs, dates, composites, **kwargs):
    ''' Get composites, such as multi-year seasonal averages or for NAO+/-

    Allows general data requests in the configured data base, e.g. ERA-Interim. The request
    can span several files, e.g. by including several vertical levels or by covering several
    years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
    
    The method internally uses metopen to locate data files. Hence, it will find data in the 
    locations given conf.datapath (in user scripts and user settings files typically 
    available as conf.datapath).

    Parameters
    ----------
    plevqs : 2-tuple or list of 2-tuples
        Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
        e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
        identifier, following the ECMWF conventions as far as appicable, e.g. ``'u'`` or ``'msl'``.
        Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
        and/or the variable name, to request all vertical levels/variables available.
    dates : list of datetime
        The minimum and maxmimum dates in this list define the requested time interval. The i
        final date will not be included in the result, i.e. for all time steps in 2016 request
        dates from 2016-01-01 00:00 to 2017-01-01 00:00.
    composites : composite_test or list of composite_test
        The compositing criteria to be used.
    
    Keyword arguments
    -----------------
    metopen arguments : all optional
        Optional arguments passed on to calls of metopen within this function.

    Returns
    -------
    dict of np.ndarray
        Composite data for the requested variable(s), vertical level(s), and composites.
    grid.gridlib
        If ``no_static=False`` meta-information about the requested data, otherwise ``None``.
    '''

    raise NotImplementedError('to be written')


# C'est le fin
