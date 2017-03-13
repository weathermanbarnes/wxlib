#!/usr/bin/env python
# -*- encoding: utf-8

from datetime import timedelta as td

from ..settings import def_context
conf = def_context('nora10')



# General settings
# ================

conf.years = range(1958,2016) # 1957 and 2016 (generally) present, but incomplete and hence skipped here.
conf.file_std = 'NORA10.%(time)s.%(plev)s.%(qf)s'
conf.file_agg = 'NORA10.%(agg)s.%(time)s.%(plev)s.%(qf)s'
conf.file_timeless = 'NORA10.%(time)s.%(name)s'
conf.file_ts = 'NORA10.%(ts)s.%(time)s.%(plev)s.%(qf)s'
conf.file_static = 'NORA10.static'
conf.timestep = td(0.125)
conf.gridsize = (400,248)
conf.datapath.insert(1, '/Data/gfi/share/Reanalyses/Met.no/NORA10')
conf.datapath.insert(1, '/Data/gfi/users/local/share')


# Variable definitions
# ====================

u = ('u', 'u', 'U component of wind', 'm s**-1')
v = ('v', 'v', 'V component of wind', 'm s**-1')
w = ('w', 'w', 'Pressure vertical velocity', 'Pa s**-1')
z = ('z', 'z', 'Geopotential', 'm**2 s**-2')
t = ('t', 't', 'Temperature', 'K')
q = ('q', 'q', 'Specific humidity', 'kg kg**-1')

msl = ('msl', 'msl', 'Sea-level pressure', 'Pa')
sp = ('sp', 'sp', 'Surface pressure', 'Pa')
u10 = ('u10', 'u10', '10 metre U wind component', 'm s**-1')
v10 = ('v10', 'v10', '10 metre V wind component', 'm s**-1')
t2m = ('t2m', 't2m', '2 metre temperature', 'K')
t0m = ('t0m', 't0m', 'Surface temperature', 'K')
slhf = ('slhf', 'slhf', 'Surface latent heat flux', 'J m**-2')
sshf = ('sshf', 'sshf', 'Surface sensible heat flux', 'J m**-2')
cp = ('cp', 'cp', 'Convective precipitation', 'mm')
lsp = ('lsp', 'lsp', 'Large-scale precipitation', 'mm')
ci = ('ci', 'ci', 'Sea-ice cover', '(0 - 1)')

oro = ('oro', None, 'Surface geopotential', 'm**2 s**-2')

# 1. Pressure levels
conf.plevs = ['100', '300', '500', '700', '850', '925', '1000', ]
conf.register_variable([u, v, w, t, z, q], conf.plevs)

# 2. Surface variables
conf.register_variable([slhf, sshf, msl, sp, u10, v10, t2m, t0m, cp, lsp, ci], ['sfc', ])

# 3. No vertical level
conf.register_variable([oro, ], [])

# that's it
