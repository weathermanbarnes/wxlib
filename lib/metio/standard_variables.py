#!/usr/bin/env python
# -*- encoding: utf-8


import numpy as np

from .datasource import variable as v


# For some variables the average is meaningless -> make histograms instead of averages
_dd_bins = [348.75, ] + list(np.arange(11.25, 360, 22.5))
_defang_bins = [17, ] + list(range(-18,18))
_defang_bins = np.array(_defang_bins)*np.pi/36.0 + np.pi/72.0


standard_variables = [
    v('u', 'U component of wind', 'm s**-1'),
    v('v', 'V component of wind', 'm s**-1'),
    v('w', 'Pressure vertical velocity', 'Pa s**-1'),
    v('ff', 'Wind speed', 'm s**-1'),
    v('dd', 'Wind direction', 'deg N', bins=_dd_bins),

    v('pv', 'Potential vorticity', 'K m**2 kg**-1 s**-1'),
    v('vo', 'Vorticity (relative)', 's**-1'),
    v('vo_curv', 'Curvature vortcity', 's**-1'),
    v('shear', 'Wind shear in natural coordinates', 's**-1'),
    v('grad_shear', 'Shear gradient in natural coordinates', 'm**-1 s**-1'),
    v('defabs', 'Total deformation', 's**-1'),
    v('defang', 'Deformation angle', 'rad', bins=_defang_bins),
    v('defanr', 'Deformation angle (natural coordinates)', 'rad', bins=_defang_bins),
    v('defstr', 'Stretching deformation in natural coordinates', 's**-1'),
    v('div', 'Divergence', 's**-1'),

    v('rsr', 'Rotation-strain ratio', '1'),
    v('ow', 'Okubo-Weiss parameter', 's**-2'),

    v('pres', 'Pressure', 'Pa'),
    v('mont', 'Montgomery potential', 'm**2 s**-2'),
    v('z', 'Geopotential', 'm**2 s**-2'),

    v('t', 'Temperature', 'K'),
    v('pt', 'Potential temperature', 'K'),
    v('eqpt', 'Equivalent potential temperature', 'K'),
    v('tdiab', 'Temperature tendency due to physics', 'K (6h)**-1'),

    v('q', 'Specific humidity', 'kg kg**-1'),
    v('lwc', 'Specific cloud liquid water content', 'kg kg**-1'),
    v('iwc', 'Specific cloud ice water content', 'kg kg**-1'),
    v('rwc', 'Specific rain water content', 'kg kg**-1'),
    v('swc', 'Specific snow water content', 'kg kg**-1'),

    v('msl', 'Mean sea level pressure', 'Pa'),
    v('sp', 'Surface pressure', 'Pa'),
    v('ci', 'Sea-ice cover', '(0 - 1)'),
    v('sst', 'Sea surface temperature', 'K'),
    v('pt2m', '2 metre potential temperature', 'K'),
    v('t2m', '2 metre temperature', 'K'),
    v('d2m', '2 metre dewpoint temperature', 'K'),
    v('u10', '10 metre U wind component', 'm s**-1'),
    v('v10', '10 metre V wind component', 'm s**-1'),
    v('ff10', '10 metre wind speed', 'm s**-1'),

    v('blh', 'Boundary layer height', 'm'),
    v('tcc', 'Total cloud cover', '1'),
    v('cape', 'Convective available potential energy', 'J kg**-1'),
    v('viewvf', 'Vertical integral of eastward water vapour flux', 'kg m**-1 s**-1'),
    v('vinwvf', 'Vertical integral of northward water vapour flux', 'kg m**-1 s**-1'),
    v('viwvd', 'Vertical integral of divergence of moisture flux', 'kg m**-2 s**-1'),
    v('vilwd', 'Vertical integral of divergence of cloud liquid water flux', 'kg m**-2 s**-1'),
    v('viiwd', 'Vertical integral of divergence of cloud frozen water flux', 'kg m**-2 s**-1'),
    v('tcw', 'Total column water', 'kg m**-2'),
    v('tcwv', 'Total column water vapour', 'kg m**-2'),
    v('tclw', 'Total column liquid water', 'kg m**-2'),
    v('tciw', 'Total column ice water', 'kg m**-2'),

    v('cp', 'Convective precipitation', 'm (6h)**-1'),
    v('csf', 'Convective snowfall', 'm of water equivalent (6h)**-1'),
    v('e', 'Evaporation', 'm of water equivalent (6h)**-1'),
    v('lsf', 'Large-scale snowfall', 'm of water equivalent (6h)**-1'),
    v('lsp', 'Large-scale precipitation', 'm (6h)**-1'),
    v('tp', 'Total precipitation', 'm (6h)**-1'),
    v('sf', 'Total snow fall', 'm of water equivalent (6h)**-1'),

    v('slhf', 'Surface latent heat flux', 'J m**-2 (6h)**-1'),
    v('sshf', 'Surface sensible heat flux', 'J m**-2 (6h)**-1'),

    v('ssr', 'Surface net solar radiation', 'J m**-2 (6h)**-1'),
    v('ssrc', 'Surface net solar radiation, clear sky', 'J m**-2 (6h)**-1'),
    v('ssrd', 'Surface downwelling solar radiation', 'J m**-2 (6h)**-1'),
    v('str', 'Surface net thermal radiation', 'J m**-2 (6h)**-1'),
    v('strc', 'Surface net thermal radiation, clear sky', 'J m**-2 (6h)**-1'),
    v('strd', 'Surface downwelling thermal radiation', 'J m**-2 (6h)**-1'),
    v('tsr', 'Top net solar radiation', 'J m**-2 (6h)**-1'),
    v('tsrc', 'Top net solar radiation, clear sky', 'J m**-2 (6h)**-1'),
    v('ttr', 'Top net thermal radiation', 'J m**-2 (6h)**-1'),
    v('ttrc', 'Top net thermal radiation, clear sky', 'J m**-2 (6h)**-1'),

    v('cold_front', 'Cold front lines', '(mixed)', lines='cold_froff'),
    v('warm_front', 'Warm front lines', '(mixed)', lines='warm_froff'),
    v('stat_front', 'Stationary front lines', '(mixed)', lines='stat_froff'),
    v('sst_front', 'SST front lines', '(mixed)', lines='sst_froff'),
    v('frovo_id', 'Frontal volume ID', '(object index)', objmask='frovo'),
    v('cold_front_freq', 'Cold front lines detection frequency', 'm of line m**-2'),
    v('warm_front_freq', 'Warm front lines detection frequency', 'm of line m**-2'),
    v('stat_front_freq', 'Stationary front lines detection frequency', 'm of line m**-2'),
    v('sst_front_freq', 'SST front lines detection frequency', 'm of line m**-2'),
    v('frovo_freq', 'Frontal volume detection frequency', '(time step)**-1'),

    v('cycmask', 'Cyclone detection mask', '(0-1)', objmask='cyc'),
    v('cyc_freq', 'Cyclone detection frequency', '(time step)**-1'),
    v('acycmask', 'Anticyclone detection mask', '(0-1)', objmask='acyc'),
    v('acyc_freq', 'Anticyclone detection frequency', '(time step)**-1'),
    v('cyc_dens', 'Cyclone detection density', '(1000 km)**-2'),
    v('cycgen_dens', 'Cyclogeneis detection density', '(1000 km)**-2'),
    v('cyclys_dens', 'Cyclolysis detection density', '(1000 km)**-2'),

    v('jetaxis', 'Jet axis lines', '(mixed)', lines='jaoff'),
    v('jetaxis_freq', 'Jet axis lines', 'm of line m**-2'),
    v('rwb_a', 'Anticyclonic wave breaking frequency', '(time step)**-1'),
    v('rwb_c', 'Cyclonic wave breaking frequency', '(time step)**-1'),
    v('block', 'Blocking frequency', '(time step)**-1'),
    v('blockint', 'Block intensity indicator', '(input) m**-1'),
]


# C'est le fin
