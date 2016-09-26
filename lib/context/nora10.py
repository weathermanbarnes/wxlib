#!/usr/bin/env python
# -*- encoding: utf-8

from datetime import timedelta as td

from ..settings import def_context
conf = def_context('nora10')



# General settings
# ================

conf.years = range(1958,2013) # 1957 and 2013 (generally) present, but incomplete and hence skipped here.
conf.file_std = 'NORA10.%(time)s.%(plev)s.%(qf)s'
conf.file_agg = 'NORA10.%(agg)s.%(time)s.%(plev)s.%(qf)s'
conf.file_timeless = 'NORA10.%(time)s.%(name)s'
conf.file_ts = 'NORA10.%(ts)s.%(time)s.%(plev)s.%(qf)s'
conf.file_static = 'NORA10.static'
conf.timestep = td(0.25)
conf.gridsize = (400,248)
conf.datapath.insert(1, '/Data/gfi/share/Reanalyses/Met.no/NORA10')
conf.datapath.insert(1, '/Data/gfi/users/local/share')


# Variable definitions
# ====================

u = ('u', 'u', 'U component of wind', 'm s**-1')
v = ('v', 'v', 'V component of wind', 'm s**-1')
w = ('w', 'w', 'Pressure vertical velocity', 'Pa s**-1')

gh = ('gh', 'gh', 'Geopotential', 'm**2 s**-2')
pt = ('pt', 'pt', 'Potential temperature', 'K')

r = ('r', 'r', 'Relative humidity', '1')

pres = ('pres', 'pres', 'Sea-level pressure', 'Pa')
u10 = ('10u', '10u', '10 metre U wind component', 'm s**-1')
v10 = ('10v', '10v', '10 metre V wind component', 'm s**-1')
t2 = ('2t', '2t', '2 metre temperature', 'K')
slhf = ('slhf', 'slhf', 'Surface latent heat flux', 'J m**-2')
sshf = ('sshf', 'sshf', 'Surface sensible heat flux', 'J m**-2')

oro = ('oro', None, 'Surface geopotential', 'm**2 s**-2')

# 1. Pressure levels
conf.plevs = ['100', '300', '500', '700', '850', '925', '1000', ]
conf.register_variable([u, v, w, pt, gh, r], conf.plevs)

# 2. No vertical level
conf.register_variable([oro, slhf, sshf, pres, u10, v10, t2], [])

# that's it
