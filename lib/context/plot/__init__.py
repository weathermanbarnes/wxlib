#!/usr/bin/env python
# -*- encoding: utf-8

import numpy as np

from ...settings import def_context, in_context, conf, cm
def_context('plot')


if in_context('erainterim'): 
	conf.plotf[None,'t','cmap'] = cm.RdBu_r
	conf.plotf[None,'pt','cmap'] = cm.RdBu_r

	conf.plotf[None,'q','cmap'] = cm.q
	conf.plotf[None,'tcw','cmap'] = cm.q
	conf.plotf[None,'tcwv','cmap'] = cm.q

	conf.plot[None,'pv','hook'] = lambda pv: pv*1e6 		# From SI units to PVU
	conf.plotf[None,'pv','hook'] = lambda pv: pv*1e6 		# -"-
	conf.plot[None, 'pv','scale'] = np.array([-2, -1, 1, 2])
	
	conf.plot[None,'q','hook'] = lambda q: q*1e3 			# From kg/kg to g/kg
	conf.plotf[None,'q','hook'] = lambda q: q*1e3 			# -"-

	conf.plot[None,'msl','hook'] = lambda msl: msl/100.0 		# From Pa to hPa
	conf.plotf[None,'msl','hook'] = lambda msl: msl/100.0 		# -"-

	conf.plot[None,'z','hook'] = lambda msl: msl/98.1		# From m2s-2 to gpdm
	conf.plotf[None,'z','hook'] = lambda msl: msl/98.1		# -"-
	
	conf.plotf[None,'oro','cmap'] = cm.gist_earth
	conf.plot[None,'oro','scale'] = range(10000,80001,10000)
	conf.plotf[None,'oro','scale'] = range(-19000,51000,2000)

if in_context('derived'):
	conf.plotf[None,'eqpt','cmap'] = cm.RdBu_r

	conf.plotf[None,'dd','scale'] = np.arange(0,36.1)*10.0
	conf.plotf[None,'dd','ticks'] = np.arange(0,8)*360.0/8.0
	conf.plotf[None,'dd','ticklabels'] = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']

	conf.plotf[None,'defabs','cmap'] = cm.defabs

	conf.plotf[None,'defang','scale'] = np.arange(-4,6)*np.pi/8 - np.pi/16 
	conf.plotf[None,'defang','colors'] = ['w', (0.625,0.625,0.85), (0.25, 0.25, 0.7), (0.375,0.375,0.6), '0.5', (0.6, 0.6, 0.375), (0.7, 0.7, 0.25), (0.85, 0.86, 0.625), 'w'] 
	conf.plotf[None,'defang','ticks'] = np.arange(-4,5)*np.pi/8.0 
	conf.plotf[None,'defang','ticklabels'] = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']
	conf.plotf[None,'defang','extend'] = 'neither'

	conf.plotf[None,'defanr','scale'] = np.arange(-4,6)*np.pi/8 - np.pi/16 
	conf.plotf[None,'defanr','colors'] = ['w', (0.625,0.625,0.85), (0.25, 0.25, 0.7), (0.375,0.375,0.6), '0.5', (0.6, 0.6, 0.375), (0.7, 0.7, 0.25), (0.85, 0.86, 0.625), 'w'] 
	conf.plotf[None,'defanr','ticks'] = np.arange(-4,5)*np.pi/8.0 
	conf.plotf[None,'defanr','ticklabels'] = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']
	conf.plotf[None,'defanr','extend'] = 'neither'

	_ticks = np.array([200,100,50,20,10])
	conf.plotf[None,'jetaxis_freq','scale'] = np.log10(1.0/_ticks)
	conf.plotf[None, 'jetaxis_freq','cmap'] = cm.defabs
	conf.plotf[None,'jetaxis_freq','ticks'] = np.log10(1.0/_ticks)
	conf.plotf[None,'jetaxis_freq','ticklabels'] = ['1/%d' % val for val in _ticks]
	conf.plotf[None,'jetaxis_freq','hook'] = lambda x: np.log10(x+1e-100)
	conf.plotf[None,'jetaxis_freq','extend'] = 'both'


