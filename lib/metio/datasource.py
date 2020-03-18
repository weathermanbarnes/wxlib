#!/usr/bin/env python
# -*- encoding: utf-8


import cftime
from datetime import timedelta as td


# 16 billion values -> 128G of memory, a bit more than 40 years of 6 hourly data at 0.5deg resolution
MAX_REQUEST_SIZE = 16.0e9 
WARN_REQUEST_SIZE = 2.0e9


timestep = None
gridsize = None
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


def metopen(filename, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
        no_static=False, quiet=False, mode='r', no_xarray=False):
    ''' Find and open files by name
    
    Uses the se.conf.datapath list to locate files in a variety of different locations. 
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

    raise NotImplementedError('Needs to be implemented per specific data source.')


def metsave(**kwargs):
    # TODO: Should may be become a generic saving function similar to xarray.to_netcdf(), 
    # but with proper group support?
    #
    # In every case: Warn about missing meta data rather than stopping with errors.
    ''' TBD '''

    raise NotImplementedError('TBD. May be data source-specific, may be not.')


# TODO: Make into a factory returning get_instantaneus tailored to the different data sets
#       (Likely this can be done through decorators?)
#       This method depends on: files_by_plevq
def get_instantaneous(plevqs, dates, force=False, **kwargs):
    ''' Get instantaneous fields

    Allows general data requests in the configured data base, e.g. ERA-Interim. The request
    can span several files, e.g. by including several vertical levels or by covering several
    years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
    
    The method internally uses metopen to locate data files. Hence, it will find data in the 
    locations given se.conf.datapath (in user scripts and user settings files typically 
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
    req = {}
    request_size = 0
    for plevq in plevqs:
        # req[plevq] contains a list of 4-tuples: 
        # (filename, list of tidx, list of dates, request_size in number of values)
        req[plevq] = list(files_by_plevq(plevq, start=start, end=end))
        request_size += sum([reqinfo[3] for reqinfo in req[plevq]])

    # Checking max length
    if request_size > MAX_REQUEST_SIZE:
        if force:
            print(f'Warning: you requested {request_size / 1024**3:.2f}G of data.')
        else:
            raise ValueError('''Cowardly refusing to fetch {request_size} values at once.
                If you are absolutely certain what you're doing, you can use force=True to override.''')

    if not force and request_size > WARN_REQUEST_SIZE:
        print(f'Warning: you requested {request_size / 1024**3:.2f}G of data.')
    
    # Remove no_static if present
    kwargs.pop('no_static', None)
    static = None

    dat = None
    for year in years:
        # Construct the slice
        fst = dt(year,1,1,0) - dt0
        lst = (dt(year+1,1,1,0) - se.conf.timestep) - dt0
        fst = int(fst.total_seconds()/dts)
        lst = int(lst.total_seconds()/dts) + 1

        # Leave out unnecessary indexes for better compatibility
        cut = slice(max(tsmin - fst, 0),min(1+tsmax - fst, lst - fst))
        datcut = slice(fst+cut.start-tsmin, fst+cut.stop-tsmin)
        
        # One or more vertical levels?
        i = 0
        if type(dat) == type(None):
            f, d, static = metopen(se.conf.file_std % {'time': year, 'plev': plevs[0], 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
            # Inject meta data for npy files
            if not f:
                static.t = np.arange(d.shape[0])*dts/3600
                static.t_parsed = [dt(year,1,1)+i*se.conf.timestep for i in range(d.shape[0])]
            if len(d.shape) == 4 and d.shape[1] > 1:
                separate_plevs = False
                s = (1+tsmax-tsmin, ) + d.shape[1:]
            elif len(d.shape) == 4: 
                separate_plevs = True
                s = (1+tsmax-tsmin, len(plevs), ) + d.shape[2:]
            else:
                separate_plevs = True
                s = (1+tsmax-tsmin, len(plevs), ) + d.shape[1:]

            dat = np.empty(s, dtype=d.dtype)
            if separate_plevs:
                if len(d.shape) < 4:
                    d = d[:,np.newaxis,::]
                dat[datcut,0:1,::] = d
            else:
                dat[datcut,::] = d

            static.t = static.t[cut]
            if type(static.t_parsed) == np.ndarray:
                static.t_parsed = static.t_parsed[cut]

            i = 1
        
        if separate_plevs:
            for plev in plevs[i:]:
                f, d, static_ = metopen(se.conf.file_std % {'time': year, 'plev': plev, 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
                if len(d.shape) < 4:
                    d = d[:,np.newaxis,::]
                dat[datcut,i:i+1,::] = d
                if i == 0:
                    static.t = np.concatenate((static.t, static_.t[cut]))
                    if type(static.t_parsed) == np.ndarray:
                        static.t_parsed = np.concatenate((static.t_parsed, static_.t_parsed[cut]))
                elif year == years[0]:
                    if not static.z_name == static_.z_name or \
                            not static.z_unit == static_.z_unit:
                        print('WARNING: concatenating vertical levels of different types!')
                        static.z_unit = u'MIXED!'
                    static.z = np.concatenate((static.z, static_.z))
                i += 1

        elif i == 0:
            f, dat[datcut,::], static_ = metopen(se.conf.file_std % {'time': year, 'plev': plevs[0], 'qf': se.conf.qf[q]}, q, cut=cut, **kwargs)
            static.t = np.concatenate((static.t, static_.t[cut]))
            if type(static.t_parsed) == np.ndarray:
                static.t_parsed = np.concatenate((static.t_parsed, static_.t_parsed[cut]))

    # Time-averaging if specified
    if tavg and len(dates) > 1:
        dat = dat.mean(axis=0)
    
    #dat = dat.squeeze()
    
    return dat, static
    raise NotImplementedError('to be written')


def get_time_average(plevqs, dates, **kwargs):
    ''' Get time-average fields

    Allows general data requests in the configured data base, e.g. ERA-Interim. The request
    can span several files, e.g. by including several vertical levels or by covering several
    years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
    
    The method internally uses metopen to locate data files. Hence, it will find data in the 
    locations given se.conf.datapath (in user scripts and user settings files typically 
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
    locations given se.conf.datapath (in user scripts and user settings files typically 
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
    locations given se.conf.datapath (in user scripts and user settings files typically 
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
