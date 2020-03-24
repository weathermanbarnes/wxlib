#!/usr/bin/env python
# -*- encoding: utf-8


import numpy as np


u = ('u', 'u', 'U component of wind', 'm s**-1')
v = ('v', 'v', 'V component of wind', 'm s**-1')
w = ('w', 'w', 'Pressure vertical velocity', 'Pa s**-1')

ff = ('ff', 'ff', 'Wind speed', 'm s**-1')
dd = ('dd', 'dd', 'Wind direction', 'deg N')

vo = ('vo', 'vo', 'Vorticity (relative)', 's**-1')
pv = ('pv', 'pv', 'Potential vorticity', 'K m**2 kg**-1 s**-1')
div = ('div', 'div', 'Divergence', 's**-1')

defabs = ('defabs', 'defabs', 'Total deformation', 's**-1')
defang = ('defang', 'defang', 'Deformation angle', 'rad N')
defanr = ('defanr', 'defanr', 'Deformation angle (natural coordinates)', 'rad')
defstr = ('defstr', 'defstr', 'Stretching deformation in natural coordinates', 's**-1')
vo_curv = ('vo_curv', 'vo_curv', 'Curvature vortcity', 's**-1')
shear = ('shear', 'shear', 'Wind shear in natural coordinates', 's**-1')
grad_shear = ('grad_shear', 'grad_shear', 'Shear gradient in natural coordinates', 'm**-1 s**-1')

rsr = ('rsr', 'rsr', 'Rotation-strain ratio', '1')
ow = ('ow', 'ow', 'Okubo-Weiss parameter', 's**-2')

pres = ('pres', 'pres', 'Pressure', 'Pa')
mont = ('mont', 'mont', 'Montgomery potential', 'm**2 s**-2')
z = ('z', 'z', 'Geopotential', 'm**2 s**-2')

t = ('t', 't', 'Temperature', 'K')
pt = ('pt', 'pt', 'Potential temperature', 'K')
eqpt = ('eqpt', 'eqpt', 'Equivalent potential temperature', 'K')
tdiab = ('tdiab', 'tdiab', 'Temperature tendency due to physics', 'K (6h)**-1')

q = ('q', 'q', 'Specific humidity', 'kg kg**-1')

msl = ('msl', 'msl', 'Mean sea level pressure', 'Pa')
sp = ('sp', 'sp', 'Surface pressure', 'Pa')
ci = ('ci', 'ci', 'Sea-ice cover', '(0 - 1)')
sst = ('sst', 'sst', 'Sea surface temperature', 'K')
pt2m = ('pt2m', 'pt2m', '2 metre potential temperature', 'K')
t2m = ('t2m', 't2m', '2 metre temperature', 'K')
tcw = ('tcw', 'tcw', 'Total column water', 'kg m**-2')
u10 = ('u10', 'u10', '10 metre U wind component', 'm s**-1')
v10 = ('v10', 'v10', '10 metre V wind component', 'm s**-1')
ff10 = ('ff10', 'ff10', '10 metre wind speed', 'm s**-1')
tcwv = ('tcwv', 'tcwv', 'Total column water vapour', 'kg m**-2')

cp = ('cp', 'cp', 'Convective precipitation', 'm (6h)**-1')
csf = ('csf', 'csf', 'Convective snowfall', 'm of water equivalent (6h)**-1')
e = ('e', 'e', 'Evaporation', 'm of water equivalent (6h)**-1')
lsf = ('lsf', 'lsf', 'Large-scale snowfall', 'm of water equivalent (6h)**-1')
lsp = ('lsp', 'lsp', 'Large-scale precipitation', 'm (6h)**-1')
tp = ('tp', 'tp', 'Total precipitation', 'm (6h)**-1')

slhf = ('slhf', 'slhf', 'Surface latent heat flux', 'J m**-2 (6h)**-1')
sshf = ('sshf', 'sshf', 'Surface sensible heat flux', 'J m**-2 (6h)**-1')

ssr = ('ssr', 'ssr', 'Surface net solar radiation', 'J m**-2 (6h)**-1')
str_ = ('str', 'str', 'Surface net thermal radiation', 'J m**-2 (6h)**-1')
tsr = ('tsr', 'tsr', 'Top net solar radiation', 'J m**-2 (6h)**-1')
ttr = ('ttr', 'ttr', 'Top net thermal radiation', 'J m**-2 (6h)**-1')

cfront = ('cold_front', 'cold_front', 'Cold front lines', '1')
cfroff = ('cold_froff', 'cold_front')
wfront = ('warm_front', 'warm_front', 'Warm front lines', '1')
wfroff = ('warm_froff', 'warm_front')
sfront = ('stat_front', 'stat_front', 'Stationary front lines', '1')
sfroff = ('stat_froff', 'stat_front')
sstfront = ('sst_front', 'sst_front', 'SST front lines', '1')
sstfroff = ('sst_froff', 'sst_front')

jetaxis = ('jetaxis', 'jetaxis', 'Jet axis lines')
jaoff = ('jaoff', 'jetaxis')

rwb_a = ('rwb_a', 'rwb_a', 'Anticyclonic wave breaking frequency', '(time step)**-1')
rwb_c = ('rwb_c', 'rwb_c', 'Cyclonic wave breaking frequency', '(time step)**-1')

blockint = ('blockint', 'blockint', 'Block intensity indicator', '(input) m**-1')
block = ('block', 'block', 'Block mask', '1')

frovo = ('frovo_id', 'frovo_id', 'Frontal volume ID', '1')
frovo_freq = ('frovo_freq', 'frovo_freq', 'Frontal volume detection frequency', '(time step)**-1')


## Some variables are more special than others

# Line detections
LINES = {
    'cold_front': 'cold_froff', 'warm_front': 'warm_froff', 'stat_front': 'stat_froff', 
    'sst_front': 'sst_froff', 
    'jetaxis': 'jaoff', 
}

# Areal detections
OBJMASK = {
    'cycmask': 'cyc',
    'frovo_id': 'frovo',
}

# Variables where the average is meaningless -> make histograms instead of averages
_dd_bins = [348.75, ] + list(np.arange(11.25, 360, 22.5))
_defang_bins = [17, ] + list(range(-18,18))
_defang_bins = np.array(_defang_bins)*np.pi/36.0 + np.pi/72.0
BINS = {
    'dd': _dd_bins,
    'defang': _defang_bins, 'defanr': _defang_bins,
}



def _average_q_name(q):
    ''' Given variable q, what is the name of its composite/time average? '''

    if q in LINES:
        qout = f'{q}_freq'
    elif q in OBJMASK:
        qout = f'{OBJMASK[q]}_freq'
    else:
        qout = q

    return qout



# C'est le fin
