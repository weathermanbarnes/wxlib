#!/usr/bin/env python
# -*- encoding: utf-8

import numpy as np

from dynlib.metio.erainterim import conf, dt, get_instantaneous, metsave
from dynlib.plotsettings import pconf, pconfl

import dynlib.diag

year = 2011
timeinterval = [dt(year,1,1,0), dt(year+1,1,1,0)]     # end of the time interval is excluced
plev = '850'

# Get wind velocity components on 850 hPa for the entire year 2011
dat, grid = get_instantaneous([(plev, 'u'), (plev,'v'), ], timeinterval)

# Calculate total deformation
defabs = dynlib.diag.def_total(dat[plev,'u'][:,0,:,:], dat[plev,'v'][:,0,:,:], grid.dx, grid.dy)

# grid does not contain vertical level information (as data from different levels could have been requested)
# -> Inject vertical level information required for saving
ogrid = grid.new_plev(plev)

# Save results as netCDF file
tosave = {
    'defabs': defabs[:,np.newaxis,:,:],      # Pressure-level data is expected to be 4-dimensional
}
metsave(tosave, ogrid, f'ei.ans.{year}.{plev}.defabs')

# the end
