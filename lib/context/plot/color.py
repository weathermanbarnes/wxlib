#!/usr/bin/env python
# -*- encoding: utf-8

from ...settings import def_context, in_context, conf, cm
def_context('plot.color')



if in_context('erainterim'):
	conf.plotf[None,'t','cmap'] = cm.RdBu_r
	conf.plotf[None,'pt','cmap'] = cm.RdBu_r

	conf.plotf[None,'q','cmap'] = cm.q
	conf.plotf[None,'tcw','cmap'] = cm.q
	conf.plotf[None,'tcwv','cmap'] = cm.q

	conf.plotf[None,'oro','cmap'] = cm.gist_earth


if in_context('derived'):
	conf.plotf[None,'defabs','cmap'] = cm.defabs
	conf.plotf[None,'defang','colors'] = ['w', (0.625,0.625,0.85), (0.25, 0.25, 0.7), (0.375,0.375,0.6), '0.5', (0.6, 0.6, 0.375), (0.7, 0.7, 0.25), (0.85, 0.86, 0.625), 'w'] 
	conf.plotf[None,'defanr','colors'] = ['w', (0.625,0.625,0.85), (0.25, 0.25, 0.7), (0.375,0.375,0.6), '0.5', (0.6, 0.6, 0.375), (0.7, 0.7, 0.25), (0.85, 0.86, 0.625), 'w'] 

	conf.plotf[None, 'jetaxis_freq','cmap'] = cm.defabs

	conf.plotf[None,'thetae','cmap'] = cm.RdBu_r


# that's it
