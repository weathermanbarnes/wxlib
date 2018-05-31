#!/usr/bin/env python
# -*- encoding: utf-8

from datetime import datetime as dt, timedelta as td

from ..settings import def_context
conf = def_context('erainterim.forecast')



# General settings
# ================

conf.epoch = dt(1979,1,1)
conf.years = range(1979,2017)
conf.file_std = 'ei.for.%(time)s.%(plev)s.%(qf)s'
conf.file_agg = 'ei.for.%(agg)s.%(time)s.%(plev)s.%(qf)s'
conf.file_timeless = 'ei.for.%(time)s.%(name)s'
conf.file_ts = 'ei.for.%(ts)s.%(time)s.%(plev)s.%(qf)s'
conf.file_static = 'ei.ans.static'
conf.timestep = td(0.25)
conf.gridsize = (361,720)
conf.datapath.insert(1, '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY')
conf.datapath.insert(1, '/Data/gfi/users/csp001/share') # for static; TODO: Move to a more general directory!


# Variable definitions
# ====================

cp = ('cp', 'cp', 'Convective precipitation', 'm (6h)**-1')
csf = ('csf', 'csf', 'Convective snowfall', 'm of water equivalent (6h)**-1')
e = ('e', 'e', 'Evaporation', 'm of water equivalent (6h)**-1')
lsf = ('lsf', 'lsf', 'Large-scale snowfall', 'm of water equivalent (6h)**-1')
lsp = ('lsp', 'lsp', 'Large-scale precipitation', 'm (6h)**-1')
slhf = ('slhf', 'slhf', 'Surface latent heat flux', 'J m**-2 (6h)**-1')
sshf = ('sshf', 'sshf', 'Surface sensible heat flux', 'J m**-2 (6h)**-1')
ssr = ('ssr', 'ssr', 'Surface net solar radiation', 'J m**-2 (6h)**-1')
str_ = ('str', 'str', 'Surface net thermal radiation', 'J m**-2 (6h)**-1')
tsr = ('tsr', 'tsr', 'Top net solar radiation', 'J m**-2 (6h)**-1')
ttr = ('ttr', 'ttr', 'Top net thermal radiation', 'J m**-2 (6h)**-1')

oro = ('oro', None, 'Surface geopotential', 'm**2 s**-2')

# 1. Surface variables
conf.sfclevs = ['sfc', ]
conf.register_variable([cp, csf, e, lsf, lsp, slhf, sshf, ssr, str_, tsr, ttr], conf.sfclevs)

# 5. No vertical level
conf.register_variable([oro, ], [])

# that's it
