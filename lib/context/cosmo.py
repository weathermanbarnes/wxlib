#!/usr/bin/env python
# -*- encoding: utf-8

from ..settings import def_context
conf = def_context('cosmo')


# Variable definitions
# ====================

u = ('U', None, 'U component of wind', 'm s**-1')
v = ('V', None, 'V component of wind', 'm s**-1')
w  = ('W', None, 'Vertical wind velocity', 'm s**-1')

t = ('T', None, 'Temperature', 'K')
ts = ('T_S', None, 'Soil surface temperature', 'Pa')
p = ('P', None, 'Pressure', 'Pa')
ps = ('PS', None, 'Surface pressure', 'Pa')
qv = ('QV', None, 'Specific humidity', 'kg kg**-1')
qvs = ('QV_S', None, 'Surface specific humidity', 'kg kg**-1')

alhfls = ('ALHFL_S', None, 'Averaged surface latent heat flux', 'W m**-2')
ashfls = ('ASHFL_S', None, 'Averaged surface sensible heat flux', 'W m**-2')

vortic_u = ('VORTIC_U', None, 'Relative vorticity, u-component in rotated grid', 's**-1')
vortic_v = ('VORTIC_V', None, 'Relative vorticity, v-component in rotated grid', 's**-1')
vortic_w = ('VORTIC_W', None, 'Relative vorticity, w-component in rotated grid', 's**-1')
pot_vortic = ('POT_VORTIC', None, 'Potential vorticity', 'K m**2 kg**-1 s**-1')
pot_t = ('POT_T', None, 'Potential temperature', 'K')

# No predefined vertical levels
conf.register_variable([u, v, w, t, ts, p, ps, qv, qvs, alhfls, ashfls, 
	vortic_u, vortic_v, vortic_w, pot_vortic, pot_t], [])


# that's it
