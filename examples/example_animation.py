#!/usr/bin/env python
# -*- encoding: utf8

# NOTE: Requires ffmpeg to be installed! Unfortunately that is neither available on skd-cyclone nor gfi-storm. 

import matplotlib as mpl
mpl.use('Agg')

from dynlib.shorthands import dt, get_instantaneous, metsave, fig, np, plt
from dynlib.settings import conf, proj

import dynlib.context.erainterim
import dynlib.context.derived

import dynlib.diag
import dynlib

from matplotlib.animation import FuncAnimation


timeinterval = [dt(2011,12,15,0), dt(2011,12,31,18)]
plev = '850'

# Get wind velocity components on 850 hPa for Ekstremværet Dagmar
u, grid = get_instantaneous('u', timeinterval, plevs=plev)
v, grid = get_instantaneous('v', timeinterval, plevs=plev)

# Calculate total deformation
defabs = dynlib.diag.def_total(u[:,0,:,:], v[:,0,:,:], grid.dx, grid.dy)

# Animate results
conf.register_variable([dynlib.context.derived.defabs, ], [plev, ])
fio = fig.setup(**conf.plotf[plev,'defabs'])

def draw_frame(tidx):
	plt.clf()
	fig.map(defabs[tidx,::], grid, q='defabs', plev=plev, 
			name=grid.t_parsed[tidx], m=proj.N_Atlantic)
	return

ani = FuncAnimation(fio, func=draw_frame, frames=range(defabs.shape[0]))
ani.save('defabs_animation.mp4', fps=4)

# the end
