#!/usr/bin/env python
# -*- encoding: utf8

from lib.shorthands import *
from lib.settings import conf
from lib import proj

import lib.context.nora10
import lib.context.plot


dat, grid = get_instantaneous('u', dt(1986,4,7,6), plevs=850)

print dat.shape
print grid.t.shape

print grid.t
print grid.t_parsed

#fig.map(dat, grid, save='nora10.png', title=None, m=proj.n_hemisphere)

lat = np.asfortranarray(grid.y.astype('f4'))
lon = np.asfortranarray(grid.x.astype('f4'))

#m = proj.n_hemisphere()
#x, y = m(grid.x, grid.y)
m = plt
x, y  = grid.x, grid.y
#m.contourf(x.flat, y.flat, dat.flat, tri=True)
#m.tricontourf(x.flat, y.flat, dat.flat)
m.contourf(dat)
#m.contourf(lon, lat, dat, latlon=True)
#m.drawparallels(range(0,90,10))
#m.drawmeridians(range(-180,180,30))
plt.colorbar()
plt.savefig('nora10.png')


#dat, grid = get_instantaneous('u', (dt(1988,1,1,0), dt(1989,4,7,6)), plevs=['300', '500', '850'])
#
#print dat.shape
#print grid.t.shape
#
#print grid.t[:5], '...', grid.t[-5:]
#print grid.t_parsed[:5], '...', grid.t_parsed[-5:]


#
