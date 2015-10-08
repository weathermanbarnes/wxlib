#!/usr/bin/env python
# -*- encoding: utf8

from dynlib.shorthands import metopen, metsave, np
from dynlib.settings import conf

import dynlib.context.erainterim
import dynlib.context.derived

import dynlib.diag


year = 2011
plev = '850'

# Get wind velocity components on 850 hPa for the entire year 2011
fu, u, grid = metopen(conf.file_std % {'qf': 'u', 'time': year, 'plev': plev}, 'u')
fv, v, grid = metopen(conf.file_std % {'qf': 'v', 'time': year, 'plev': plev}, 'v')

# Calculate total deformation
defabs = dynlib.diag.def_total(u[:,0,:,:], v[:,0,:,:], grid.dx, grid.dy)

# Save results as netCDF file
metsave(defabs[:,np.newaxis,:,:], grid, q='defabs', plev=plev)

# the end
