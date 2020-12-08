#!/usr/bin/env python
# -*- encoding: utf-8


from .datasource import *
from ..interpol import vert_by_coord
from ..settings import settings_obj, default_conf


conf = settings_obj({
    'q': {},
    'qf': {}, 
    'qstd': {},
    'q_units': {},
    'q_long': {},
    'q_bins': {},
    'datapath': ['.',
        '/Data/gfi/share/era5/ml',
        '/Data/gfi/users/local/share',
        '/Data/gfi/users/csp001/share', 
    ], 
    'opath': '.',
    'oformat': 'nc',
    'plotpath': '.',
    'plotformat': 'png',
    'staticfile': 'ei.ans.static',
    'epoch': dt(1900,1,1,0),
    'calendar': 'standard',
    'timestep': td(0.125),
    'gridsize': (361,720),
    'local_timezone': default_conf.local_timezone,
}, [])
dt = cftime.DatetimeGregorian
# TODO: Register variables and move the LINES, OBJMASK, BINS definitions back to the conf object


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

average_q_name  = average_q_name_factory()

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
            
            if q in Hq and (plev in Hplev or plev == '__all_plev__' or plev == '__all__'):
                self.filetype = 'H'
                if plev in Hplev:
                    self.nlevs = 1
                else:
                    self.nlevs = len(Hplev)
            elif q in Zq and (plev in Zplev or plev == '__sel_plev__' or plev == '__all__'):
                self.filetype = 'Z'
                if plev in Zplev:
                    self.nlevs = 1
                else:
                    self.nlevs = len(Zplev)
            elif q in Nq:
                self.filetype = 'N'
                self.nlevs = 1
            elif q in Pq:
                self.filetype = 'P'
                if plev == '__all__' or plev == '__all_mlev__':
                    self.nlevs = n_mlevs         # All available model levels
                elif plev == '__all_plev__': 
                    self.nlevs = len(Hplev)
                elif plev == '__sel_plev__':
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
            size = (1, self.nlevs, ) + gridsize

            self.cur += timestep
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
    
    fo, oro = metopen(staticfile, 'oro', verbose=verbose, no_dtype_conversion=no_dtype_conversion, 
            quiet=quiet, no_static=True)
    static = grid_by_static(fo)
    static.oro = oro[::]
    static.t_epoch = dt(1900,1,1,0)
    static.t_unit = 'hours since 1900-01-01 00:00:00'
    static.t_interval_unit = 3600
    fo.close()

    return static


# Derive data source-specific version of metopen
metopen = metopen_factory(get_static)
metsave, metsave_composite = metsave_factory(metopen)


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
            elif plev in ['__all__', '__all__mlevs__']:
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


# Derive data source-specific versions of the remaining data getter functions
get_instantaneous = get_instantaneous_factory(files_by_plevq, metopen, get_from_file)
get_time_average = get_time_average_factory(files_by_plevq, get_from_file, get_static)
get_aggregate = get_aggregate_factory(files_by_plevq, get_from_file, get_static)
get_composite = get_composite_factory(files_by_plevq, get_from_file, get_static)


# C'est le fin
