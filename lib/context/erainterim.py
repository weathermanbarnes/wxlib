#!/usr/bin/env python
# -*- encoding: utf-8

from datetime import timedelta as td

from ..settings import def_context, conf
def_context('erainterim')



# General settings
# ================

conf.years = range(1979,2015)
conf.file_std = 'ei.ans.%(time)s.%(plev)s.%(qf)s'
conf.file_agg = 'ei.ans.%(agg)s.%(time)s.%(plev)s.%(qf)s'
conf.file_timeless = 'ei.ans.%(time)s.%(name)s'
conf.file_ts = 'ei.ans.%(ts)s.%(time)s.%(plev)s.%(qf)s'
conf.file_static = 'ei.ans.static'
conf.timestep = td(0.25)
conf.gridsize = (361,720)
conf.datapath.insert(1, '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY')
conf.datapath.insert(1, '/Data/gfi/users/csp001/share') # for static; TODO: Move to a more general directory!


# Variable definitions
# ====================

u = ('u', 'u', 'U component of wind', 'm s**-1')
v = ('v', 'v', 'V component of wind', 'm s**-1')
w = ('w', 'w', 'Pressure vertical velocity', 'Pa s**-1')

vo = ('vo', 'zeta', 'Vorticity (relative)', 's**-1')
pv = ('pv', 'pv', 'Potential vorticity', 'K m**2 kg**-1 s**-1')

pres = ('pres', 'p', 'Pressure', 'Pa')
mont = ('mont', 'm', 'Montgomery potential', 'm**2 s**-2')
z = ('z', 'Z', 'Geopotential', 'm**2 s**-2')
t = ('t', 'T', 'Temperature', 'K')
pt = ('pt', 'th', 'Potential temperature', 'K')

q = ('q', 'q', 'Specific humidity', 'kg kg**-1')

msl = ('msl', 'msl', 'Mean sea level pressure', 'Pa')
sp = ('sp', 'ps', 'Surface pressure', 'Pa')
ci = ('ci', 'SIC', 'Sea-ice cover', '(0 - 1)')
sst = ('sst', 'SST', 'Sea surface temperature', 'K')
t2m = ('t2m', 'T2', '2 metre temperature', 'K')
tcw = ('tcw', 'tw', 'Total column water', 'kg m**-2')
u10 = ('u10', 'u10', '10 metre U wind component', 'm s**-1')
v10 = ('v10', 'v10', '10 metre V wind component', 'm s**-1')
tcwv = ('tcwv', 'wv', 'Total column water vapour', 'kg m**-2')

oro = ('oro', None, 'Surface geopotential', 'm**2 s**-2')

# 1. Pressure levels
conf.plevs = ['100', '200', '300', '400', '500', '550', '600', '650', '700', '750', '800', '850', '900', '950', '1000', ]
conf.register_variable([u, v, w, pv, z, t, q], conf.plevs)

# 2. Potential temperature levels
conf.ptlevs = ['pt300', 'pt315', 'pt330', 'pt350', ]
conf.register_variable([u, v, vo, pv, pres, mont, z], conf.ptlevs)

# 3. PV-levels
conf.pvlevs = ['pv2000', ]
conf.register_variable([u, v, pt, z], conf.pvlevs)

# 4. Surface variables
conf.sfclevs = ['sfc', ]
conf.register_variable([msl, sp, ci, sst, t2m, tcw, u10, v10, tcwv], conf.sfclevs)

# 5. No vertical level
conf.register_variable([oro, ], [])

# that's it
