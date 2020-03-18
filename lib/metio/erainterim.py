#!/usr/bin/env python
# -*- encoding: utf-8


from datasource import *


timestep = td(0.25)
gridsize = (361,720)


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
            size= len(tidxs) * gridsize[0] * gridsize[1]

            self.cur = dt(self.cur.year+1, 1, 1, 0)
            return filename, tidxs, dates, size

        else:
            raise StopIteration


# C'est le fin
