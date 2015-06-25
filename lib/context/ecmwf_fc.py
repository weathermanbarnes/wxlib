#!/usr/bin/env python
# -*- encoding: utf-8

import os

from ..settings import def_context, conf
def_context('ecmwf_fc')



# General settings
# ================

_PATH = '/Data/gfi/metno/ecmwf'
conf.times = [t for t in os.listdir(_PATH) if os.path.isdir(os.path.join(_PATH,t))]
conf.file_std = '%(time)d/ec.for.%(time)d.%(plev)s.%(q)s'
conf.datapath.insert(1, _PATH)

# Variable definitions
# ====================

u = ('u', 'u', 'U component of wind', 'm s**-1')
v = ('v', 'v', 'V component of wind', 'm s**-1')
w = ('w', 'w', 'Vertical velocity', 'Pa s**-1')

pv = ('pv', 'PV', 'Potential vorticity', 'K m**2 kg**-1 s**-1')

z = ('z', 'Z', 'Geopotential', 'm**2 s**-2')
t = ('t', 'T', 'Temperature', 'K')
pt = ('pt', 'th', 'Potential temperature', 'K')

q = ('q', 'q', 'Specific humidity', 'kg kg**-1')

msl = ('msl', 'msl', 'Mean sea level pressure', 'Pa')
ci = ('ci', 'SIC', 'Sea-ice cover', '(0 - 1)')
sst = ('sst', 'SST', 'Sea surface temperature', 'K')
t2m = ('t2m', 'T2', '2 metre temperature', 'K')
d2m = ('d2m', 'TD2', '2 metre dewpoint temperature', 'K')
u10 = ('u10', 'u10', '10 metre U wind component', 'm s**-1')
v10 = ('v10', 'v10', '10 metre V wind component', 'm s**-1')
tcw = ('tcw', 'tw', 'Total column water', 'kg m**-2')
tcwv = ('tcwv', 'tv', 'Total column water vapour', 'kg m**-2')
cape = ('cape', 'CAPE', 'Convective available potential energy', 'J kg**-1')
cin = ('cin', 'CIN', 'Convective inhibition', 'J kg**-1')
cp = ('cp', 'cprec', 'Convective precipitation', 'm')
lsp = ('lsp', 'lprec', 'Large-scale precipitation', 'm')
hcc = ('hcc', 'hcl', 'High cloud cover', '(0 - 1)')
mcc = ('mcc', 'mcl', 'Medium cloud cover', '(0 - 1)')
lcc = ('lcc', 'lcl', 'Low cloud cover', '(0 - 1)')
tcc = ('tcc', 'tcl', 'Total cloud cover', '(0 - 1)')
slhf = ('slhf', 'lhf', 'Surface latent heat flux', 'J m**-2')
sshf = ('sshf', 'shf', 'Surface sensible heat flux', 'J m**-2')


# 1. Pressure levels
conf.plevs = ['300', '500', '700', '850', '925', ]
conf.register_variable([u, v, w, z, t, q], conf.plevs)

# 2. Potential temperature levels
conf.ptlevs = ['pt300', 'pt315', 'pt330', 'pt350', ]
conf.register_variable([u, v, pv], conf.ptlevs)

# 3. PV-levels
conf.pvlevs = ['pv2000', ]
conf.register_variable([u, v, pt], conf.pvlevs)

# 4. Surface variables
conf.sfclevs = ['sfc', ]
conf.register_variable([msl, ci, sst, t2m, d2m, u10, v10, tcw, tcwv, cape, cin, 
	cp, lsp, hcc, mcc, lcc, tcc, slhf, sshf], conf.sfclevs)


# that's it
