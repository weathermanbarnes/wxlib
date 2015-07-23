#!/usr/bin/env python
# -*- encoding: utf-8

from datetime import timedelta as td

from ..settings import def_context, conf
def_context('ncep1')



# General settings
# ================

conf.years = range(1948,2012) # 2012 files present but do not cover entire year, hence excluded here.
conf.file_std = '%(qf)s.%(time)s'
conf.file_agg = 'nc.ans.%(agg)s.%(time)s.%(qf)s'
conf.file_timeless = 'nc.ans.%(time)s.%(name)s'
conf.file_ts = 'nc.ans.%(ts)s.%(time)s.%(plev)s.%(qf)s'
conf.file_static = 'nc.ans.static'
conf.timestep = td(0.25)
conf.gridsize = (73,144)
conf.datapath.insert(1, '/Data/gfi/share/Reanalyses/NCEP1/6_HOURLY/PRESSURE_LEVELS')
conf.datapath.insert(1, '/Data/gfi/share/Reanalyses/NCEP1/6_HOURLY/SURFACE')
conf.datapath.insert(1, '/Data/gfi/users/local/share')


# Variable definitions
# ====================

uwnd = ('uwnd', 'uwnd', 'U component of wind', 'm s**-1')
vwnd = ('vwnd', 'vwnd', 'V component of wind', 'm s**-1')
omega = ('omega', 'omega', 'Pressure vertical velocity', 'Pa s**-1')

air = ('air', 'air', 'Temperature', 'K')
shum = ('shum', 'shum', 'Specific humidity', 'kg kg**-1')

slp = ('slp', 'slp', 'Mean sea level pressure', 'Pa')
pres = ('pres', 'pres.sfc', 'Surface pressure', 'Pa')

pr_wtr = ('pr_wtr', 'pr_wtr.eatm', 'Precipitable water', 'kg m**-2')


oro = ('oro', None, 'Surface geopotential', 'm**2 s**-2')


# 1. Pressure levels
conf.plevs = ['10', '20', '30', '50', '70', '100', '150', '200', '300', '400', '500', '600', '700', '850', '925', '1000', ]
conf.register_variable([uwnd, vwnd, omega, air, shum], conf.plevs)

# 2. Surface variables
conf.sfclevs = ['sfc', ]
conf.register_variable([slp, pres, pr_wtr], conf.sfclevs)

# 3. No vertical level
conf.register_variable([oro, ], [])

# that's it
