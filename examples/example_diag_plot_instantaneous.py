#!/usr/bin/env python
# -*- encoding: utf8

from dynlib.shorthands import dt, get_instantaneous, metsave, fig, np
from dynlib.settings import conf, proj

import dynlib.context.erainterim
import dynlib.context.derived

import dynlib.diag


timeinterval = [dt(2011,12,25,0), dt(2011,12,26,12)]
plev = '850'

# Get wind velocity components on 850 hPa for Ekstremværet Dagmar
u, grid = get_instantaneous('u', timeinterval, plevs=plev)
v, grid = get_instantaneous('v', timeinterval, plevs=plev)

# Calculate total deformation
defabs = dynlib.diag.def_total(u[:,0,:,:], v[:,0,:,:], grid.dx, grid.dy)

# Save results as netCDF file
metsave(defabs[:,np.newaxis,:,:], grid, q='defabs', plev=plev)

# Plot results
conf.register_variable([dynlib.context.derived.defabs, ], [plev, ])
for tidx in range(len(grid.t)):
	fig.map(defabs[tidx,::], grid, q='defabs', plev=plev, 
			name=grid.t_parsed[tidx], m=proj.N_Atlantic)

# the end
