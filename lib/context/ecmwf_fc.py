#!/usr/bin/env python
# -*- encoding: utf-8

import os
import os.path

from ..settings import def_context
conf = def_context('ecmwf_fc')



# General settings
# ================

_PATH = '/Data/gfi/metno/ecmwf'
if os.path.exists(_PATH):
	conf.times = [t for t in os.listdir(_PATH) if os.path.isdir(os.path.join(_PATH,t))]
	conf.datapath.insert(1, _PATH)
conf.file_std = '%(time)d/ec.for.%(time)d.%(plev)s.%(q)s'
conf.file_ens = '%(time)d/ec.ens.%(time)d.%(plev)s.%(q)s'
conf.file_ensctrl = '%(time)d/ec.ensctrl.%(time)d.%(plev)s.%(q)s'

# Variable definitions
# ====================

u = ('u', 'u', 'U component of wind', 'm s**-1')
v = ('v', 'v', 'V component of wind', 'm s**-1')
w = ('w', 'w', 'Vertical velocity', 'Pa s**-1')

pv = ('pv', 'pv', 'Potential vorticity', 'K m**2 kg**-1 s**-1')

z = ('z', 'z', 'Geopotential', 'm**2 s**-2')
t = ('t', 't', 'Temperature', 'K')
pt = ('pt', 'pt', 'Potential temperature', 'K')

q = ('q', 'q', 'Specific humidity', 'kg kg**-1')

msl = ('msl', 'msl', 'Mean sea level pressure', 'Pa')
ic = ('siconc', 'ic', 'Sea-ice concentration', '(0 - 1)')
sst = ('sst', 'sst', 'Sea surface temperature', 'K')
t2m = ('t2m', 't2m', '2 metre temperature', 'K')
d2m = ('d2m', 'd2m', '2 metre dewpoint temperature', 'K')
u10 = ('u10', 'u10', '10 metre U wind component', 'm s**-1')
v10 = ('v10', 'v10', '10 metre V wind component', 'm s**-1')
tcw = ('tcw', 'tcw', 'Total column water', 'kg m**-2')
tcwv = ('tcwv', 'tcwv', 'Total column water vapour', 'kg m**-2')
cape = ('cape', 'cape', 'Convective available potential energy', 'J kg**-1')
cin = ('cin', 'cin', 'Convective inhibition', 'J kg**-1')
cp = ('cp', 'cp', 'Convective precipitation', 'm')
lsp = ('lsp', 'lsp', 'Large-scale precipitation', 'm')
hcc = ('hcc', 'hcc', 'High cloud cover', '(0 - 1)')
mcc = ('mcc', 'mcc', 'Medium cloud cover', '(0 - 1)')
lcc = ('lcc', 'lcc', 'Low cloud cover', '(0 - 1)')
tcc = ('tcc', 'tcc', 'Total cloud cover', '(0 - 1)')
slhf = ('slhf', 'slhf', 'Surface latent heat flux', 'J m**-2')
sshf = ('sshf', 'sshf', 'Surface sensible heat flux', 'J m**-2')


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
conf.register_variable([msl, ic, sst, t2m, d2m, u10, v10, tcw, tcwv, cape, cin, 
	cp, lsp, hcc, mcc, lcc, tcc, slhf, sshf], conf.sfclevs)


# that's it
