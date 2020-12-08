#!/usr/bin/env python
# -*- encoding: utf-8


from .datasource import *
from .standard_variables import LINES, OBJMASK, BINS
from ..settings import settings_obj, default_conf


conf = settings_obj({
    'q': {},
    'qf': {}, 
    'qstd': {},
    'q_units': {},
    'q_long': {},
    'q_bins': {},
    'datapath': ['.',
        '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY', 
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
    'timestep': td(0.25),
    'gridsize': (361,720),
    'local_timezone': default_conf.local_timezone,
}, [])
dt = cftime.DatetimeGregorian
# TODO: Register variables and move the LINES, OBJMASK, BINS definitions back to the conf object


average_q_name = average_q_name_factory()

_files_by_plevq = files_by_plevq
class files_by_plevq(_files_by_plevq):
    def __init__(self, plevq, start=None, end=None):
        plev, q = plevq
        if plev == '__all__':
            raise ValueError('ERA-Interim does not support requests for all vertical levels.')
        if q == '__all__':
            raise ValueError('ERA-Interim does not support requests for all variables.')

        if not start:
            start = dt(1979,1,1,0)
        if not end:
            end = dt(2019,1,1,0)
        
        return super().__init__(plevq, start, end)

    def __next__(self):
        if self.cur < self.end:
            filename = f'ei.ans.{self.cur.year}.{self.plev}.{self.q}'

            yearlen = int((dt(self.cur.year+1, 1, 1, 0) - dt(self.cur.year, 1, 1, 0)).total_seconds() / timestep.total_seconds())
            tidxs_all = range(yearlen)
            dates_all = [dt(self.cur.year, 1, 1, 0) + td(0.25)*i for i in tidxs_all]

            tidxs = [tidx for tidx in tidxs_all if dates_all[tidx] >= self.start and dates_all[tidx] < self.end]
            dates = [dates_all[tidx] for tidx in tidxs_all if dates_all[tidx] >= self.start and dates_all[tidx] < self.end]
            size = (len(tidxs), 1) + gridsize

            self.cur = dt(self.cur.year+1, 1, 1, 0)
            return filename, tidxs, dates, size

        else:
            raise StopIteration


def get_static(verbose=False, no_dtype_conversion=False, quiet=False):
    ''' Get standard meta-information for ERA-Interim

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

    if plev not in filename:
        raise ValueError(f'Filename `{filename}` does not match the requested vertical level `{plev}.`')
    
    # Pass everything on, but the file object itself.
    stuff = metopen(filename, q, **kwargs)
    if len(stuff) == 2:
        return stuff[1]
    else:
        return stuff[1:]


# Derive data source-specific versions of the remaining data getter functions
get_instantaneous = get_instantaneous_factory(files_by_plevq, metopen, get_from_file)
get_time_average = get_time_average_factory(files_by_plevq, get_from_file, get_static)
get_aggregate = get_aggregate_factory(files_by_plevq, get_from_file, get_static)
get_composite = get_composite_factory(files_by_plevq, get_from_file, get_static)


# C'est le fin
