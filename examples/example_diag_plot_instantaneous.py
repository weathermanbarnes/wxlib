#!/usr/bin/env python
# -*- encoding: utf-8

import numpy as np

from dynlib.metio.erainterim import conf, dt, get_instantaneous, metsave
from dynlib.plotsettings import pconf, pconfl

import dynlib.diag
import dynlib.proj as proj
import dynlib.figures as fig

timeinterval = [dt(2011,12,25,0), dt(2011,12,26,12)]
plev = '850'

# Get wind velocity components on 850 hPa for Ekstremværet Dagmar
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
metsave(tosave, ogrid, f'ei.ans.Dagmar.{plev}.defabs')

# Plot results
for tidx in range(len(grid.t)):
    fig.map(defabs[tidx,::], grid, q='defabs', plev=plev, 
            name=grid.t_parsed[tidx], m=proj.N_Atlantic, show=False,
            save=f'Dagmar_defabs_{grid.t_parsed[tidx].strftime("%Y%m%d_%H")}.pdf')
    fig.plt.close()

# the end
