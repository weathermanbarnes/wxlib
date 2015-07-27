#!/usr/bin/env python
# -*- encoding: utf8

from lib.shorthands import *
from lib.settings import conf

import lib.context.erainterim
import lib.context.plot


dat, grid = get_instantaneous('u', dt(1986,4,7,6), plevs=850)

print dat.shape
print grid.t.shape

print grid.t
print grid.t_parsed

fig.map(dat, grid, save='erai.png', title=None)


dat, grid = get_instantaneous('u', (dt(1988,1,1,0), dt(1989,4,7,6)), plevs=['300', '500', '850'])

print dat.shape
print grid.t.shape

print grid.t[:5], '...', grid.t[-5:]
print grid.t_parsed[:5], '...', grid.t_parsed[-5:]


#
