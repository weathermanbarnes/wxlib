#!/usr/bin/env python
# -*- encoding: utf-8

import numpy as np
from ..settings import def_context
conf = def_context('derived', parent='active')




# Variable definitions
# ====================

ff = ('ff', 'ff', 'Wind speed', 'm s**-1')
dd = ('dd', 'dd', 'Wind direction', '(0 - 360)')

vo = ('vo', 'vo', 'Vorticity (relative)', 's**-1')
div = ('div', 'div', 'Divergence', 's**-1')
defabs = ('defabs', 'defabs', 'Total deformation', 's**-1')
_bins_defang = [17,]
_bins_defang.extend(range(-18,18))
_bins_defang = np.array(_bins_defang)*np.pi/36.0 + np.pi/72.0
defang = ('defang', 'defang', 'Deformation angle', 'rad', _bins_defang)
defanr = ('defanr', 'defanr', 'Deformation angle (natural coordinates)', 'rad', _bins_defang)

rsr = ('rsr', 'rsr', 'Rotation-strain ratio', '1')
ow = ('ow', 'ow', 'Okubo-Weiss parameter', 's**-2')

pt = ('pt', 'pt', 'Potential temperature', 'K')
eqpt = ('eqpt', 'eqpt', 'Equivalent potential temperature', 'K')

cfront = ('cold_front', 'cold_front', 'Cold front lines', '1')
cfroff = ('cold_froff', 'cold_front')
wfront = ('warm_front', 'warm_front', 'Warm front lines', '1')
wfroff = ('warm_froff', 'warm_front')
sfront = ('stat_front', 'stat_front', 'Stationary front lines', '1')
sfroff = ('stat_froff', 'stat_front')
sstfront = ('sst_front', 'sst_front', 'SST front lines', '1')
sstfroff = ('sst_froff', 'sst_front')
cfront_freq = ('cold_front_freq', 'cold_front_freq', 'Cold front detection frequency', 'm m**-2')
wfront_freq = ('warm_front_freq', 'warm_front_freq', 'Warm front detection frequency', 'm m**-2')
sfront_freq = ('stat_front_freq', 'stat_front_freq', 'Stationary front detection frequency', 'm m**-2')
sstfront_freq = ('sst_front_freq', 'sst_front_freq', 'SST front detection frequency', 'm m**-2')

vorl = ('vorl', 'vorl', 'Vorticity lines', '1')
vloff = ('vloff', 'vorl')
convl = ('convl', 'convl', 'Convergence lines', '1')
cloff = ('cloff', 'convl')
defl = ('defl', 'defl', 'Deformation lines', '1')
dloff = ('dloff', 'defl')

shear = ('shear', 'shear', 'Wind shear in natural coordinates', 's**-1')
grad_shear = ('grad_shear', 'grad_shear', 'Shear gradient in natural coordinates', 'm**-1 s**-1')
jetaxis = ('jetaxis', 'jetaxis', 'Jet axis lines')
jaoff = ('jaoff', 'jetaxis')
jetaxis_freq = ('jetaxis_freq', 'jetaxis_freq', 'Jet axis detection frequency', '(time step)**-1')

rwb_a = ('rwb_a', 'rwb_a', 'Anticyclonic wave breaking frequency', '(time step)**-1')
rwb_c = ('rwb_c', 'rwb_c', 'Cyclonic wave breaking frequency', '(time step)**-1')

blockint = ('blockint', 'blockint', 'Block intensity indicator', '(input) m**-1')
block = ('block', 'block', 'Block mask', '1')
block_freq = ('block_freq', 'block_freq', 'Block detection frequency', '(time step)**-1')

frovo = ('frovo_id', 'frovo_id', 'Frontal volume ID', '1')
frovo_freq = ('frovo_freq', 'frovo_freq', 'Frontal volume detection frequency', '(time step)**-1')

cyc = ('cycmask', 'cycmask', 'Cyclone object ID', '1')
cyc_freq = ('cyc_freq', 'cyc_freq', 'Cyclone detection frequency', '(time step)**-1')
cyc_dens = ('cyc_dens', 'cyc_dens', 'Cyclone detection density', '(time step)**-1 (1000 km)**-1')
cycgen_dens = ('cycgen_dens', 'cycgen_dens', 'Cyclogeneis detection density', '(time step)**-1 (1000 km)**-1')
cyclys_dens = ('cyclys_dens', 'cyclys_dens', 'Cyclolysis detection density', '(time step)**-1 (1000 km)**-1')

# The vertical level of the following variables is fixed and does not depend on user application
conf.register_variable([sstfront, sstfroff, sstfront_freq, cyc, cyc_freq, cyc_dens, cycgen_dens, cyclys_dens], ['sfc',])

# The vertical levels on which the following variables are available depends 
# on the application, and must hence be defined in the user settings.
conf.register_variable([ff, dd, vo, div, defabs, defang, defanr, rsr, ow, 
	pt, eqpt, frovo, frovo_freq,
	cfront, cfroff, cfront_freq, wfront, wfroff, wfront_freq, sfront, sfroff, sfront_freq,
	vorl, vloff, convl, cloff, defl, dloff, 
	grad_shear, jetaxis, jaoff, jetaxis_freq, 
	rwb_a, rwb_c, blockint, block, block_freq], [])


# that's it
