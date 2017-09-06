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

cfront = ('cfront', 'cfront', 'Cold front lines')
cfroff = ('cfroff', 'cfront')
wfront = ('wfront', 'wfront', 'Warm front lines')
wfroff = ('wfroff', 'wfront')
sfront = ('sfront', 'sfront', 'Stationary front lines')
sfroff = ('sfroff', 'sfront')

vorl = ('vorl', 'vorl', 'Vorticity lines')
vloff = ('vloff', 'vorl')
convl = ('convl', 'convl', 'Convergence lines')
cloff = ('cloff', 'convl')
defl = ('defl', 'defl', 'Deformation lines')
dloff = ('dloff', 'defl')

grad_shear = ('grad_shear', 'grad_shear', 'Shear gradient in natural coordinates', 'm**-1 s**-1')
jetaxis = ('jetaxis', 'jetaxis', 'Jet axis lines')
jaoff = ('jaoff', 'jetaxis')
jetaxis_freq = ('jetaxis_freq', 'jetaxis_freq', 'Jet axis detection frequency', '(time step)**-1')

rwb_a = ('rwb_a', 'rwb_a', 'Anticyclonic wave breaking frequency', '(time step)**-1')
rwb_c = ('rwb_c', 'rwb_c', 'Cyclonic wave breaking frequency', '(time step)**-1')

blockint = ('blockint', 'blockint', 'Block intensity indicator', '(input) m**-1')
block = ('block', 'block', 'Block mask', '1')


#conf.new_variable('div', 'defabs_tend_div')
#conf.new_variable('defabs_tend_beta', 'defabs_tend_beta')
#conf.new_variable('defabs_tend_prescor', 'defabs_tend_prescor')
#conf.new_variable('defabs_tend_tilt', 'defabs_tend_tilt')
#conf.new_variable('defabs_tend_adv', 'defabs_tend_adv')
#conf.new_variable('defabs_tend_adv3d', 'defabs_tend_adv3d')

#conf.new_variable(None, 'defabs_tend_lagrange')
#conf.new_variable(None, 'defabs_tend_euler')

#conf.new_variable('defang_tend_beta', 'defang_tend_beta')
#conf.new_variable('defang_tend_prescor', 'defang_tend_prescor')
#conf.new_variable('defang_tend_tilt', 'defang_tend_tilt')
#conf.new_variable('defang_tend_adv', 'defang_tend_adv')
#conf.new_variable('defang_tend_adv3d', 'defang_tend_adv3d')



# The vertical levels on which any of these variables are available depends 
# on the application, and must hence be defined in the user settings.
conf.register_variable([ff, dd, vo, div, defabs, defang, defanr, rsr, ow, 
	pt, eqpt, 
	cfront, cfroff, wfront, wfroff, sfront, sfroff, 
	vorl, vloff, convl, cloff, defl, dloff, 
	grad_shear, jetaxis, jaoff, jetaxis_freq, 
	rwb_a, rwb_c, blockint, block], [])


# that's it
