#!/usr/bin/env python
# -*- encoding: utf-8


from .datasource import *
from ..interpol import vert_by_coord
from ..settings import settings_obj, default_conf

import scipy.interpolate as interp


dt = cftime.DatetimeGregorian
conf = settings_obj({
    'q': {},
    'qf': {}, 
    'q_std': {},
    'q_avg': {},
    'q_units': {},
    'q_long': {},
    'q_bins': {},
    'q_lines': {},
    'q_obj': {},
    'datapath': ['.',
        '/Data/gfi/share/era5/ml',
        '/Data/gfi/users/local/share',
        '/Data/gfi/users/csp001/share', 
    ], 
    'opath': '.',
    'oformat': 'nc',
    'staticfile': 'ei.ans.static',
    'epoch': dt(1900,1,1,0),
    'calendar': 'standard',
    'timestep': td(0.125),
    'gridsize': (361,720),
    'local_timezone': default_conf.local_timezone,
}, [])

# Renaming for convenience
v = variable
variables = [
    v('U', 'U component of wind', 'm s**-1', q_std='u'),
    v('V', 'V component of wind', 'm s**-1', q_std='v'),
    v('OMEGA', 'Pressure vertical velocity', 'Pa s**-1', q_std='w'),

    v('PV', 'Potential vorticity', 'K m**2 kg**-1 s**-1', q_std='pv'),
    v('Z', 'Geopotential', 'm**2 s**-2', q_std='z'),

    v('T', 'Temperature', 'K', q_std='t'),
    v('PT', 'Potential temperature', 'K', q_std='pt'),
    v('EQPT', 'Equivalent potential temperature', 'K', q_std='eqpt'),

    v('Q', 'Specific humidity', 'kg kg**-1', q_std='q'),
    v('LWC', 'Specific cloud liquid water content', 'kg kg**-1', q_std='lwc'),
    v('IWC', 'Specific cloud ice water content', 'kg kg**-1', q_std='iwc'),
    v('RWC', 'Specific rain water content', 'kg kg**-1', q_std='rwc'),
    v('SWC', 'Specific snow water content', 'kg kg**-1', q_std='swc'),

    v('MSL', 'Mean sea level pressure', 'Pa', q_std='msl'),
    v('PS', 'Surface pressure', 'Pa', q_std='sp'),
    v('CI', 'Sea-ice cover', '(0 - 1)', q_std='ci'),
    v('SSTK', 'Sea surface skin temperature', 'K', q_std='sst'),
    v('T2M', '2 metre temperature', 'K', q_std='t2m'),
    v('D2M', '2 metre dewpoint temperature', 'K', q_std='d2m'),
    v('U10M', '10 metre U wind component', 'm s**-1', q_std='u10'),
    v('V10M', '10 metre V wind component', 'm s**-1', q_std='v10'),
    v('WG10', 'Wind gust at 10 metres', 'm s**-1'),

    v('BLH', 'Boundary layer height', 'm', q_std='blh'),
    v('TCC', 'Total cloud cover', '1', q_std='tcc'),
    v('VIWVD', 'Vertical integral of divergence of moisture flux', 'kg m**-2 s**-1', q_std='viwvd'),
    v('TCW', 'Total column water', 'kg m**-2', q_std='tcw'),
    v('TCWV', 'Total column water vapour', 'kg m**-2', q_std='tcwv'),
    v('TCLW', 'Total column liquid water', 'kg m**-2', q_std='tclw'),
    v('TCIW', 'Total column ice water', 'kg m**-2', q_std='tciw'),

    v('CP', 'Convective precipitation', 'm (6h)**-1', q_std='cp'),
    v('E', 'Evaporation', 'm of water equivalent (6h)**-1', q_std='e'),
    v('LSP', 'Large-scale precipitation', 'm (6h)**-1', q_std='lsp'),
    v('SF', 'Total snow fall', 'm of water equivalent (6h)**-1', q_std='sf'),

    v('SLHF', 'Surface latent heat flux', 'W m**-2', q_std='slhf'),
    v('SSHF', 'Surface sensible heat flux', 'W m**-2', q_std='sshf'),

    v('SI', 'Solar insolation', 'W m**-2'),
    v('SSR', 'Surface net solar radiation', 'W m**-2', q_std='ssr'),
    v('SSRC', 'Surface net solar radiation, clear sky', 'W m**-2', q_std='ssrc'),
    v('SSRD', 'Surface downwelling solar radiation', 'W m**-2', q_std='ssrd'),
    v('STR', 'Surface net thermal radiation', 'W m**-2', q_std='str'),
    v('STRC', 'Surface net thermal radiation, clear sky', 'W m**-2', q_std='strc'),
    v('STRD', 'Surface downwelling thermal radiation', 'W m**-2', q_std='strd'),
    v('TSR', 'Top net solar radiation', 'W m**-2', q_std='tsr'),
    v('TSRC', 'Top net solar radiation, clear sky', 'W m**-2', q_std='tsrc'),
    v('TTR', 'Top net thermal radiation', 'W m**-2', q_std='ttr'),
    v('TTRC', 'Top net thermal radiation, clear sky', 'W m**-2', q_std='ttrc'),
]
conf.register_variable(variables)


filetypes = ['B', 'H', 'N', 'P', 'Z']
# Variables present in B-files (surface variables --- disregarded in favour of the N files)
Bq = [ 'PS', 'MSL', 'TCC', 'U10M', 'V10M', 'SSTK', 'CI', 'D2M', 'T2M', 'TCW', 'TCWV', 'VIWVD', 'E', 
        'MN2T', 'MX2T', 'SI', 'TTR', 'TTRC', 'WG10', 'LSP', 'CP', 'SF', 'SSHF', 'SLHF', 'BLH', ]
# Variables present in Z-files (geopotential on 37 pressure levels)
Hq = ['Z', ]
Hplev = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125',
         '150', '175', '200', '225', '250', '300', '350', '400', '450', '500',
         '550', '600', '650', '700', '750', '775', '800', '825', '850', '875',
         '900', '925', '950', '975', '1000', ]
# Variables present in N-files (surface variables, extended set)
Nq = [ 'PS', 'MSL', 'TCC', 'U10M', 'V10M', 'SSTK', 'CI', 'D2M', 'T2M', 'TCW', 'TCWV', 'VIWVD', 'E', 
        'MN2T', 'MX2T', 'SI', 'TTR', 'TTRC', 'WG10', 'LSP', 'CP', 'SF', 'SSHF', 'SLHF', 'BLH', 
        'SSR', 'STR', 'SSRD', 'STRD', 'TSR', 'TSRC', 'SSRC', 'STRC', 'TCLW', 'TCIW', 'CAPE', 'FAL', ]
# Variables present in P-files (primary 3D variables on model levels and few selected surface variables)
Pq = [ 'LSP', 'CP', 'SF', 'SSHF', 'SLHF', 'BLH', 'PS',
        'Q',  'T', 'OMEGA', 'LWC', 'IWC', 'RWC', 'SWC', 'U', 'V', ]
n_mlevs = 98
# Variables present in Z-files (selected variables on 11 selected pressure levels)
Zq = ['U', 'V', 'Q', 'T', 'Z', ]
Zplev = ['100', '200', '250', '300', '400', '500', '600', '700', '800', '850', '900']


OBJMASK = {}
LINES = {}
BINS = {}

_files_by_plevq = files_by_plevq
class files_by_plevq(_files_by_plevq):
    def __init__(self, plevq=None, filetype='', start=None, end=None):
        if filetype:
            if filetype not in filetypes:
                raise ValueError(f'Unknown filetype `{filetype}`, must be one of {filetypes}.')

            self.filetype = filetype
            plevq = ('__all__', '__all__')

            if filetype == 'H':
                self.nlevs = len(Hplev)
            elif filetype == 'P':
                self.nlevs = n_mlevs
            elif filetype == 'Z':
                self.nlevs = len(Zplev)
            else:
                self.nlevs = 1

        elif plevq:
            plev, q = plevq
            if q == '__all__':
                raise ValueError('ERA-5 does not support requests for all variables.'
                        'If you want to request all variables in a given file type, '
                        'supply the filetype argument instead of plevq.')
            
            if q in Hq and (plev in Hplev or plev == '__all_plevs__' or plev == '__all__'):
                self.filetype = 'H'
                if plev in Hplev:
                    self.nlevs = 1
                else:
                    self.nlevs = len(Hplev)
            elif q in Zq and (plev in Zplev or plev == '__sel_plevs__'):
                self.filetype = 'Z'
                if plev in Zplev:
                    self.nlevs = 1
                else:
                    self.nlevs = len(Zplev)
            elif q in Bq:                       # Temporary fix for N files being incomplete.
                self.filetype = 'B'
                self.nlevs = 1
            elif q in Nq:
                self.filetype = 'N'
                self.nlevs = 1
            elif q in Pq:
                self.filetype = 'P'
                if plev == '__all__' or plev == '__all_mlevs__':
                    self.nlevs = n_mlevs         # All available model levels
                elif plev == '__all_plevs__': 
                    self.nlevs = len(Hplev)
                elif plev == '__sel_plevs__':
                    self.nlevs = len(Zplev)
                else:                       # Everything else will be intepreted as a single level specification
                    self.nlevs = 1
            else:
                raise ValueError(f'Unknown variable for requested vertical level (`{q}` on `{plev}`).')
        
        else:
            raise ValueError('Either plevq or filetype must be provided as arguments.')

        if not start:
            start = dt(1980,1,1,0)
        if not end:
            end = dt(2020,1,1,0)
        
        return super().__init__(plevq, start, end)

    def __next__(self):
        if self.cur < self.end:
            filename = f'{self.cur.year}/{self.cur.month:02d}/{self.filetype}{self.cur.year}{self.cur.month:02d}{self.cur.day:02d}_{self.cur.hour:02d}'
            dates = [self.cur, ]
            tidxs = [0, ]
            size = (1, self.nlevs, ) + conf.gridsize

            self.cur += conf.timestep
            return filename, tidxs, dates, size

        else:
            raise StopIteration


def get_static(verbose=False, no_dtype_conversion=False, quiet=False):
    ''' Get standard meta-information for ERA-5

    Parameters
    ----------
    quiet : bool
        *Optional*, default ``False``. Suppress default output from metopen.
    verbose : bool
        *Optional*, default ``False``. Print debug information on where the static
        file is being looked for.
    no_dtype_conversion : bool
        *Optional*, default ``False``. By default, ``metopen`` uncompresses data in the 
        file and converts all data to float64. This behaviour can be suppressed by 
        setting this option to ``True``.

    Returns
    -------
    gridlib.grid
        Some meta information about the data, like the grid information.
    '''
    
    fo, oro = metopen(conf.staticfile, 'oro', verbose=verbose, no_dtype_conversion=no_dtype_conversion, 
            quiet=quiet, no_static=True)
    static = grid_by_static(fo)
    static.oro = oro[::]
    static.t_epoch = dt(1900,1,1,0)
    static.t_unit = 'hours since 1900-01-01 00:00:00'
    static.t_interval_unit = 3600
    fo.close()

    return static


# Derive data source-specific version of metopen
metopen = metopen_factory(get_static, conf)
metsave, metsave_composite = metsave_factory(metopen, conf)


_get_from_file = get_from_file
def get_from_file(filename, plev, q, **kwargs):
    __doc__ = _get_from_file.__doc__
    
    # Which kind of file are we dealing with?
    filetype = filename[8:9]
    stuff = metopen(filename, 'time', **kwargs)
    f = stuff[0]

    # Model level files in general require vertical interpolation
    if filetype in ['P','S']:
        # There is 98 model levels for the 3D variables, but a and b are defined on the original 137 model levels.
        # TODO: I assume the 98 chosen levels are the 98 lowest ones; verify whether this is correct!
        # TODO: I assume the metadata for PS is wrong, and this is actually hPa rather than Pa
        ps = f['PS'][::]
        a = f['hyam'][-n_mlevs:]
        b = f['hybm'][-n_mlevs:]
        p = a[np.newaxis,:,np.newaxis,np.newaxis] + b[np.newaxis,:,np.newaxis,np.newaxis] * ps[:,np.newaxis,:,:]*100
        if len(stuff) == 3:
            stuff[2].hybrid_p = p 

        if q == '__all__':
            dat = {}
            for q in globals()[f'{filetype}q']:
                datr = f[q][::]
                if len(datr.shape) == 3:
                    dat[q] = datr
                elif plev in ['__all__', '__all_mlevs__']:
                    dat[q] = datr
                else:
                    dat[q] = vert_by_coord(datr, p, np.array([int(plev)*100, ]) )

        else: 
            datr = f[q][::]
            if len(datr.shape) == 3:
                dat = datr
            elif plev in ['__all__', '__all_mlevs__']:
                dat = datr
            else:
                dat = vert_by_coord(datr, p, np.array([int(plev)*100, ]) ) 

    # Extract the correct vertical level if only one is requested, otherwise return all available plevs
    elif filetype in ['H', 'Z']:
        avail_plevs = globals()[f'{filetype}plev']
        if plev in avail_plevs:
            pidx = avail_plevs.index(plev)
            pslc = slice(pidx,pidx+1)
        elif plev in ['__all__', '__all_plevs__', '__sel_plevs__']:
            pslc = slice(None)
        else:
            raise ValueError('plev must be either a single pressure level, or one of the special strings '
                    f'`__all__`, `__all_plevs__`, `__sel_plevs__`. Got `{plev}` instead.')

        if q == '__all__':
            dat = {}
            for q in globals()[f'{filetype}q']:
                dat[q] = f[q][:,pslc,:,:]

        else:
            dat = f[q][:,pslc,:,:]

    # Surface data: More-or-less straightforward read and pass on, everything except the file object
    else:
        if q == '__all__':
            dat = {}
            for q in globals()[f'{filetype}q']:
                dat[q] = f[q][::]

        else:
            dat = f[q][::]

    f.close()
    
    if len(stuff) == 3:
        return dat, stuff[2]
    else:
        return dat

get_normalized_from_file = get_normalized_from_file_factory(get_from_file, conf)

# Derive data source-specific versions of the remaining data getter functions
get_instantaneous = get_instantaneous_factory(files_by_plevq, metopen, get_from_file, get_static, conf)
get_time_average = get_time_average_factory(files_by_plevq, get_normalized_from_file, get_static, conf)
get_aggregate = get_aggregate_factory(files_by_plevq, get_normalized_from_file, get_static, conf)
get_composite = get_composite_factory(files_by_plevq, get_normalized_from_file, get_static, conf)


def get_hor_interpolation_functions(dates, plevs, q):
    ''' Get values of the variable q at given dates, interpolated to given plevs, ys, xs

    The given plevs, ys and xs must all have the same 4-dimensional shape. Dates must be
    one-dimensional, corresponding in length to the first dimension of plevs, ys and xs.
    The returned array will have the same shape as plevs, ys and xs.
    
    Parameters
    ----------
    dates : list/np.ndarray of datetime with dimensions (N,)
        The dates for which data is to be interpolated.
    plevs : np.ndarray with dimensions (z,)
        Pressure levels to be interpolated to.
    q : str
        A variable name identifier, following the ECMWF conventions as far as applicable,
        e.g. ``'u'`` or ``'msl'``.

    Returns
    -------
    dict (date => dict (plev => callable))
        A nested dictionary of interpolation functions, with the requested dates and plevs as keys.
    '''

    # Assuming the chunking of data into files is consistent across plevs, this is a natural 
    # rhythm to iterate through both test and composite data
    #
    # (the below code does not require data to be chunked consistently, 
    #  but the code will much more efficient if it is)
    if type(dates) == list:
        dates = np.array(dates)
    if type(plevs) == list:
        plevs = np.array(plevs)

    start, end = dates.min(), dates.max() + td(0,1) # the last requested date should be included here.
    if len(plevs) > 1:
        req_plevq = ('__all_mlevs__', q) 
    else:
        req_plevq = (plevs[0], q)

    req = list(files_by_plevq(req_plevq, start=start, end=end))
    
    # Initialise dictionary to hold the results
    ifuncs = {}

    # Iterate through the data set in chunks natural to the data source
    for filename, tidxs, dates_, shape in req:
        # Skip chunks from which no time step is required
        load_chunk = False
        for date in dates_:
            if date in dates:
                load_chunk = True

        if not load_chunk:
            continue
        
        dat_, grid = get_normalized_from_file(filename, *req_plevq)
        if len(plevs) > 1:
            dat__ = dat_
            dat_ = np.empty(dat_.shape[:1] + (len(plevs),) + dat_.shape[2:])
            for tidx in range(dat_.shape[0]):
                dat_[tidx,...] = vert_by_coord(dat__[tidx,...], grid.hybrid_p, plevs*100)
        
        for tidx_in, date in enumerate(dates_):
            tidxs_out = np.argwhere(dates == date)
            if tidxs_out.size == 0:
                continue
            
            for pidx, plev in enumerate(plevs):
                if len(dat_.shape) == 3:
                    if pidx > 0:
                        raise ValueError('Expected several levels in dat_ array, but only got 3D-array')
                    dat__ = dat_[tidx_in,:,:]
                else:
                    dat__ = dat_[tidx_in,pidx,:,:]
                    
                if np.isnan(dat__).sum() > 0:
                    zonalavg = np.nanmean(dat__, axis=1)
                    dat__ -= zonalavg[:,np.newaxis]
                    dat__, conv = utils.fill_nan(dat__[np.newaxis,:,:])
                    dat__ = dat__[0,:,:]
                    dat__ += zonalavg[:,np.newaxis]

                if not date in ifuncs:
                    ifuncs[date] = {}
                ifuncs[date][plev] = interp.RectBivariateSpline(grid.y[:,0], grid.x[0,:], dat__, kx=1, ky=1)

    return ifuncs



def get_at_position(dates, plevs, ys, xs, q):
    ''' Get values of the variable q at given dates, interpolated to given plevs, ys, xs

    The given plevs, ys and xs must all have the same 4-dimensional shape. Dates must be
    one-dimensional, corresponding in length to the first dimension of plevs, ys and xs.
    The returned array will have the same shape as plevs, ys and xs.
    
    Parameters
    ----------
    dates : list/np.ndarray of datetime with dimensions (N,)
        The dates for which data is to be interpolated.
    plevs : np.ndarray with dimensions (z,)
        Pressure levels to be interpolated to.
    ys : np.ndarray with dimensions (N,*grid)
        Meridional distances/latitudes to be interpolated to.
    xs : np.ndarray with dimensions (N,*grid)
        Zonal distances or longitudes to be interpolated to.
    q : str
        A variable name identifier, following the ECMWF conventions as far as applicable,
        e.g. ``'u'`` or ``'msl'``.

    Returns
    -------
    np.ndarray with dimensions (N,z,*grid)
        Interpolated values of variable q at the givens dates/positions.
    '''

    # TODO: This is an awful amount of duplicated code with the above
    # TODO: ... and actually also with the generic variant in dynlib.metio.datasource

    # Consistency checks
    if not xs.shape == ys.shape:
        raise ValueError('xs and ys must have the same shape.')
        
    if not xs.shape[0] == len(dates):
        raise ValueError('Number of dates must be equal length of first dimension of xs and ys.')
    
    # Assuming the chunking of data into files is consistent across plevs, this is a natural 
    # rhythm to iterate through both test and composite data
    #
    # (the below code does not require data to be chunked consistently, 
    #  but the code will much more efficient if it is)
    if type(dates) == list:
        dates = np.array(dates)

    start, end = dates.min(), dates.max() + td(0,1) # the last requested date should be included here.
    if len(plevs) > 1:
        req_plevq = ('__all_mlevs__', q) 
    else:
        req_plevq = (plevs[0], q)

    req = list(files_by_plevq(req_plevq, start=start, end=end))
    
    # Allocate array to hold the results
    dat = np.empty(xs.shape[:1] + (len(plevs),) + (xs.shape[1:]))

    # Iterate through the data set in chunks natural to the data source
    for filename, tidxs, dates_, shape in req:
        # Skip chunks from which no time step is required
        load_chunk = False
        for date in dates_:
            if date in dates:
                load_chunk = True

        if not load_chunk:
            continue
        
        dat_, grid = get_normalized_from_file(filename, *req_plevq)
        if len(plevs) > 1:
            dat__ = dat_
            dat_ = np.empty(dat_.shape[:1] + (len(plevs),) + dat_.shape[2:])
            for tidx in range(dat_.shape[0]):
                dat_[tidx,...] = vert_by_coord(dat__[tidx,...], grid.hybrid_p, plevs*100)
        
        for tidx_in, date in enumerate(dates_):
            tidxs_out = np.argwhere(dates == date)
            if tidxs_out.size == 0:
                continue
            
            for pidx in range(dat_.shape[1]):
                if len(dat_.shape) == 3:
                    if pidx > 0:
                        raise ValueError('Expected several levels in dat_ array, but only got 3D-array')
                    dat__ = dat_[tidx_in,:,:]
                else:
                    dat__ = dat_[tidx_in,pidx,:,:]
                    
                if np.isnan(dat__).sum() > 0:
                    zonalavg = np.nanmean(dat__, axis=1)
                    dat__ -= zonalavg[:,np.newaxis]
                    dat__, conv = utils.fill_nan(dat__[np.newaxis,:,:])
                    dat__ = dat__[0,:,:]
                    dat__ += zonalavg[:,np.newaxis]

                ifunc = interp.RectBivariateSpline(grid.y[:,0], grid.x[0,:], dat__, kx=1, ky=1)
                for tidx_out in tidxs_out:
                    tidx_out = tidx_out[0]
                    dat[tidx_out,pidx,...] = ifunc(ys[tidx_out,...], xs[tidx_out,...], grid=False)

    return dat


# C'est le fin
