#!/usr/bin/env python
# -*- encoding: utf-8

from ..settings import def_context, conf
def_context('wrf')


# Variable definitions
# ====================

uu = ('UU', None, 'U component of wind', 'm s**-1')
vv = ('VV', None, 'V component of wind', 'm s**-1')
w  = ('W', None, 'W component of wind', 'm s**-1')
ww = ('WW', None, 'Pressure vertical velocity', 'Pa s**-1')

tt = ('TT', None, 'Temperature', 'K')
ght = ('GHT', None, 'Geopotential', 'm')
pres = ('PRES', None, 'Pressure', 'Pa')
mu = ('MU', None, 'Perturbation dry air mass in column', 'Pa')
mub = ('MUB', None, 'Base state dry air mass in column', 'Pa')
psfc = ('PSFC', None, 'Surface pressure', 'Pa')

qvapor = ('QVAPOR', None, 'Water vapor mixing ratio', 'kg kg**-1')
qcloud = ('QCLOUD', None, 'Cloud water mixing ratio', 'kg kg**-1')
qrain = ('QRAIN', None, 'Rain water mixing ratio', 'kg kg**-1')
qice = ('QICE', None, 'Ice mixing ratio', 'kg kg**-1')
qsnow = ('QSNOW', None, 'Snow mixing ratio', 'kg kg**-1')
qgraup = ('QGRAUP', None, 'Graupel mixing ratio', 'kg kg**-1')

rh = ('RH', None, 'Relative humidity', '%')
cldfra = ('CLDFRA', None, 'Cloud fraction', '(0 - 1)')
rainc = ('RAINC', None, 'Accumulated total cumulus precipitation', 'mm')
rainnc = ('RAINNC', None, 'Accumulated total grid-scale precipitation', 'mm')

h_diabatic = ('H_DIABATIC', None, 'Latent heating due to microphysics', 'K s**-1')
rthcuten = ('RTHCUTEN', None, 'Coupled potential temperature tendency due to cumulus scheme', 'Pa K s**-1')
rthblten = ('RTHBLTEN', None, 'Coupled potential temperature tendency due to boundary layer parameterization', 'Pa K s**-1')

pblh = ('PBLH', None, 'Boundary layer height', 'm')

hfx = ('HFX', None, 'Upward heat flux at the surface', 'W m**-2')
qfx = ('QFX', None, 'Upward moisture flux at the surface', 'kg m**-2 s**-1')
lh = ('LH', None, 'Latent heat flux at the surface', 'W m**-2')


# No predefined vertical levels
conf.register_variable([uu, vv, w, ww, tt, ght, pres, mu, mub, psfc, 
	qvapor, qcloud, qrain, qice, qsnow, qgraup, rh, cldfra, rainc, rainnc, 
	h_diabatic, rthcuten, rthblten, pblh, hfx, qfx, lh], [])


# that's it
