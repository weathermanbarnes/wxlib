#!/usr/bin/env python
# -*- encoding: utf-8


from .datasource import *
from .standard_variables import standard_variables
from ..settings import settings_obj, default_conf


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
        '/Data/gfi/share/Reanalyses/NORA10', 
        '/Data/gfi/users/local/share',
    ], 
    'opath': '.',
    'oformat': 'nc',
    'staticfile': 'NORA10.197001.sfc.msl.nc',
    'epoch': dt(1900,1,1,0),
    'calendar': 'standard',
    'timestep': td(0.125),
    'gridsize': (400,248),
    'local_timezone': default_conf.local_timezone,
}, [])
# NORA10 is following our standard naming convention
conf.register_variable(standard_variables)


_files_by_plevq = files_by_plevq
class files_by_plevq(_files_by_plevq):
    def __init__(self, plevq, start=None, end=None):
        plev, q = plevq
        if plev == '__all__':
            raise ValueError('NORA10 does not support requests for all vertical levels.')
        if q == '__all__':
            raise ValueError('NORA10 does not support requests for all variables.')

        if not start:
            start = dt(1979,1,1,0)
        if not end:
            end = dt(2019,1,1,0)
        
        return super().__init__(plevq, start, end)

    def __next__(self):
        if self.cur < self.end:
            filename = f'NORA10.{self.cur.year}{self.cur.month:02d}.{self.plev}.{conf.qf.get(self.q, self.q)}'
            
            nxt = dt(self.cur.year + (self.cur.month // 12), self.cur.month % 12 + 1, 1, 0)
            monlen = int((
                nxt - dt(self.cur.year, self.cur.month, 1, 0)
            ).total_seconds() / conf.timestep.total_seconds())
            tidxs_all = range(monlen)
            dates_all = [dt(self.cur.year, self.cur.month, 1, 0) + conf.timestep*i for i in tidxs_all]

            tidxs = [tidx for tidx in tidxs_all if dates_all[tidx] >= self.start and dates_all[tidx] < self.end]
            dates = [dates_all[tidx] for tidx in tidxs_all if dates_all[tidx] >= self.start and dates_all[tidx] < self.end]

            if self.plev == 'sfc':
                size = (len(tidxs),) + conf.gridsize
            else:
                size = (len(tidxs), 1) + conf.gridsize

            self.cur = nxt
            return filename, tidxs, dates, size

        else:
            raise StopIteration


def get_static(verbose=False, no_dtype_conversion=False, quiet=False):
    ''' Get standard meta-information for NORA10

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
    
    fo, static = metopen(conf.staticfile, verbose=verbose, no_dtype_conversion=no_dtype_conversion, quiet=quiet)
    static.oro = np.zeros(conf.gridsize)

    static = static.new_time([])
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
get_time_average = get_time_average_factory(files_by_plevq, get_normalized_from_file, get_static, conf)
get_aggregate = get_aggregate_factory(files_by_plevq, get_normalized_from_file, get_static, conf)
get_composite = get_composite_factory(files_by_plevq, get_normalized_from_file, get_static, conf)

# Does not work for NORA10 as lat/lon are not monotonously increasing in the NORA10 grid directions
#get_at_position, get_hor_interpolation_functions = get_at_position_factory(files_by_plevq, get_normalized_from_file)


# C'est le fin
