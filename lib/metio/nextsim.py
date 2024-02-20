#!/usr/bin/env python
# -*- encoding: utf-8


from .datasource import *
from .standard_variables import standard_variables
from ..settings import settings_obj, default_conf
from ..utils import scale


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
    'q_feature_dists': {},
    'q_obj': {},
    'datapath': ['.',
        '/Data/gfi/spengler/csp001/seaice_tendencies',
    ], 
    'opath': '.',
    'oformat': 'nc',
    'staticfile': 'nextsim.static.t_grid.nc',
    'epoch': dt(1900,1,1,0),
    'calendar': 'standard',
    'timestep': td(0.25),
    'gridsize': (603,528),
    'local_timezone': default_conf.local_timezone,
}, [])


_files_by_plevq = files_by_plevq
class files_by_plevq(_files_by_plevq):
    def __init__(self, plevq, start=None, end=None):
        plev, q = plevq
        if not start:
            start = dt(1995,1,1,0)
        if not end:
            end = dt(2019,1,1,0)
        
        return super().__init__(plevq, start, end)

    def __next__(self):
        if self.cur < self.end:
            filename = f'nextsim.{self.cur.year}{self.cur.month:02d}.sfc.seaice_tendencies'
            
            nxt = dt(self.cur.year + (self.cur.month // 12), self.cur.month % 12 + 1, 1, 0)
            monlen = int((
                nxt - dt(self.cur.year, self.cur.month, 1, 0)
            ).total_seconds() / conf.timestep.total_seconds())
            tidxs_all = range(monlen)
            dates_all = [dt(self.cur.year, self.cur.month, 1, 0) + conf.timestep*i for i in tidxs_all]

            tidxs = [tidx for tidx in tidxs_all if dates_all[tidx] >= self.start and dates_all[tidx] < self.end]
            dates = [dates_all[tidx] for tidx in tidxs_all if dates_all[tidx] >= self.start and dates_all[tidx] < self.end]

            size = (len(tidxs), ) + conf.gridsize

            self.cur = nxt
            return filename, tidxs, dates, size

        else:
            raise StopIteration


def get_static(verbose=False, no_dtype_conversion=False, quiet=False):
    ''' Get standard meta-information for ERA5

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
    
    fo, lsm, static = metopen(conf.staticfile, 'landmask', verbose=verbose, no_dtype_conversion=no_dtype_conversion, 
            quiet=quiet)
    
    # Derive grid information from the variables in the static file
    static.x = scale(fo['longitude'])
    static.y = scale(fo['latitude'])
    static.dx = scale(fo['dx']) * 2
    static.dy = scale(fo['dy']) * 2

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

    if plev not in filename:
        raise ValueError(f'Filename `{filename}` does not match the requested vertical level `{plev}.`')
    
    # Pass everything on, but the file object itself.
    stuff = metopen(filename, q, **kwargs)

    if len(stuff) == 2:
        return stuff[1]
    else:
        return stuff[1:]

get_normalized_from_file = get_normalized_from_file_factory(get_from_file, conf)

# Derive data source-specific versions of the remaining data getter functions
get_instantaneous = get_instantaneous_factory(files_by_plevq, metopen, get_from_file, get_static, conf)
get_instantaneous_normalized = get_instantaneous_factory(files_by_plevq, metopen, get_normalized_from_file, get_static, conf, normalize_output=True)
get_time_average = get_time_average_factory(files_by_plevq, get_normalized_from_file, get_static, conf)
get_aggregate = get_aggregate_factory(files_by_plevq, get_normalized_from_file, get_static, conf)
get_composite = get_composite_factory(files_by_plevq, get_normalized_from_file, get_static, conf)

get_at_position, get_hor_interpolation_functions, get_hif_nans = get_at_position_factory(files_by_plevq, get_normalized_from_file)


# C'est le fin
