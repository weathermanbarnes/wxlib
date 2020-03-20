#!/usr/bin/env python
# -*- encoding: utf-8


from .datasource import *


timestep = td(0.25)
gridsize = (361,720)
staticfile = 'ei.ans.static'

conf.datapath.insert(1, '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY')
conf.datapath.insert(1, '/Data/gfi/users/csp001/share') # for static; TODO: Move to a more general directory!



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


_get_from_file = get_from_file
def get_from_file(filename, plev, q, **kwargs):
    __doc__ = _get_from_file.__doc__

    if plev not in filename:
        raise ValueError(f'Filename `{filename}` does not match the requested vertical level `{plev}.`')

    return metopen(filename, q, **kwargs)


# Derive data source-specific versions of the remaining data getter functions
get_instantaneous = get_instantaneous_factory(files_by_plevq, get_from_file, get_static)

# C'est le fin
