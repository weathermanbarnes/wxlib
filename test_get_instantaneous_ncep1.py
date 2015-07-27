#!/usr/bin/env python
# -*- encoding: utf8

from lib.shorthands import *
from lib.settings import conf

import lib.context.ncep1
import lib.context.plot


dat, grid = get_instantaneous('uwnd', dt(1986,4,7,6))

print dat.shape
print grid.t.shape

print grid.t
print grid.t_parsed

fig.map(dat[2], grid, save='ncep1.png', title=None)


dat, grid = get_instantaneous('uwnd', (dt(1948,1,1,0), dt(1949,4,7,6)))

print dat.shape
print grid.t.shape

print grid.t[:5], '...', grid.t[-5:]
print grid.t_parsed[:5], '...', grid.t_parsed[-5:]



#
