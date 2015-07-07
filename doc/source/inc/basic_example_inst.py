#!/usr/bin/env python
# -*- encoding: utf8

from dynlib.shorthands import dt, get_instantaneous, metsave, fig
from dynlib.settings import conf, m

import dynlib.context.erainterim
import dynlib.context.derived

import dynlib.diag


timeinterval = [dt(2011,12,25,0), dt(2011,12,26,12)]
plev = '850'

# Get wind velocity components on 850 hPa for Ekstremværet Dagmar
u, grid = get_instantaneous('u', timeinterval, plevs=plev)
v, grid = get_instantaneous('v', timeinterval, plevs=plev)

# Calculate total deformation
defabs = dynlib.diag.def_total(u, v, grid.dx, grid.dy)

# Save results as netCDF file
# TODO: Currently time and plev are not part of/taken from grid object
metsave(defabs, grid, q='defabs')

# Plot results
for tidx in range(len(grid.t)):
	fig.map(defabs[tidx,::], grid, q='defabs', plev=plev, 
			name=grid.t_parsed[tidx], proj=m.N_Atlantic)

# the end
