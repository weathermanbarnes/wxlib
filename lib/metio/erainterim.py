#!/usr/bin/env python
# -*- encoding: utf-8


from datasource import *


timestep = td(0.25)
gridsize = (361,720)


_files_by_plev_q = files_by_plev_q
class files_by_plev_q(_files_by_plev_q):
    def __init__(self, plev, q, start=None, end=None):
        if not start:
            start = dt(1979,1,1,0)
        if not end:
            end = dt(2019,1,1,0)
        
        return super().__init__(plev, q, start, end)

    def __next__(self):
        if self.cur < self.end:
            filename = f'ei.ans.{self.cur.year}.{self.plev}.{self.q}'
            self.cur = dt(self.cur.year+1,1,1,0)
            return filename
        else:
            raise StopIteration


# C'est le fin
