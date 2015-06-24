#!/usr/bin/env python
# -*- encoding: utf-8

from ..settings import def_context, conf
def_context('derived')



# Variable definitions
# ====================

ff = ('ff', 'ff', 'Wind speed', 'm s**-1')
dd = ('dd', 'dd', 'Wind direction', '(0 - 360)')

vo = ('vo', 'zeta', 'Vorticity (relative)', 's**-1')
div = ('div', 'div', 'Divergence', 's**-1')
defabs = ('defabs', 'defabs', 'Total deformation', 's**-1')
defang = ('defang', 'defang', 'Deformation angle', 'rad')
defanr = ('defanr', 'defanr', 'Deformation angle (natural coordinates)', 'rad')

rsr = ('rsr', 'rsr', 'Rotation-strain ratio', '1')
ow = ('ow', 'ow', 'Okubo-Weiss parameter', 's**-2')

eqpt = ('eqpt', 'the', 'Equivalent potential temperature', 'K')

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

jetaxis = ('jetaxis', 'jetaxis', 'Jet axis lines')
jaoff = ('jaoff', 'jetaxis')
jetaxis_freq = ('jetaxis_freq', 'jetaxis_freq', 'Jet axis detection frequency', '(time step)**-1')

rwb_a = ('rwb_a', 'rwb_a', 'Anticyclonic wave breaking frequency', '(time step)**-1')
rwb_c = ('rwb_c', 'rwb_c', 'Cyclonic wave breaking frequency', '(time step)**-1')

blockint = ('blockint', 'blockint', 'Block intensity indicator')
block = ('block', 'block', 'Block mask', '1')


# The vertical levels on which any of these variables are available depends 
# on the application, and must hence be defined in the user settings.
conf.register_variable([ff, dd, vo, div, defabs, defang, defanr, rsr, ow, 
	eqpt, 
	cfront, cfroff, wfront, wfroff, sfront, sfroff, 
	vorl, vloff, convl, cloff, defl, dloff, 
	jetaxis, jaoff, rwb_a, rwb_c, blockint, block], [])


# that's it
