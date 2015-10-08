#!/usr/bin/env python
# -*- encoding: utf8

from dynlib.shorthands import dt, get_instantaneous, metsave, metsave_lines, fig, np, dynfor
from dynlib.settings import conf, proj

import dynlib.context.erainterim
import dynlib.context.derived

import dynlib
import dynlib.detect
import dynlib.utils


timeinterval = [dt(2011,12,25,0), dt(2011,12,26,12)]
plev = 'pv2000'


# Get wind velocity components on 2 PVU for Ekstremværet Dagmar
u, grid = get_instantaneous('u', timeinterval, plevs=plev)
v, grid = get_instantaneous('v', timeinterval, plevs=plev)

# Detect jets
dynfor.config.grid_cyclic_ew = grid.cyclic_ew
dynfor.config.jetint_thres = 5.5e-9
jetaxis, jaoff = dynlib.detect.jetaxis(10000, 5000, u[:,0,:,:], v[:,0,:,:], grid.dx, grid.dy)

# Save lines as netCDF file
metsave_lines(jetaxis, jaoff, grid, q='jetaxis', qoff='jaoff', plev=plev)

# Create a grid mask of line locations and save this mask
s = grid.dx.shape
jetaxis_mask = dynlib.utils.mask_lines(s[1], s[0], jetaxis, jaoff)
conf.register_variable([dynlib.context.derived.jetaxis, dynlib.context.derived.jetaxis_freq, ], ['pv2000', ])
metsave(jetaxis_mask[:,np.newaxis,:,:], grid, q='jetaxis_freq', plev=plev)

# Plot results and overlay jet axis lines
ff = np.sqrt(u*u + v*v)
conf.register_variable([dynlib.context.derived.ff, ], [plev, ])
for tidx in range(len(grid.t)):
	fig.setup(**conf.plotf[plev,'ff'])
	overlays = [
		fig.map_overlay_lines(jetaxis[tidx,:,:], jaoff[tidx,:], grid, q='jetaxis', plev='pv2000', linecolor='k')
	]
	fig.map(ff[tidx,0,:,:], grid, q='ff', plev=plev, overlays=overlays, save='test_%02d.png' % tidx,
			name=grid.t_parsed[tidx], m=proj.N_Atlantic)

# the end
