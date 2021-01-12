#!/usr/bin/env python
# -*- encoding: utf-8

import numpy as np
import matplotlib as mpl

from . import pconf, pconfl
from .. import cm


###############################################
# Generic plot defaults 
pconf [None,None,'coastcolor'] = 'k'
pconf [None,None,'coastwidth'] = 0.5
pconf [None,None,'gridcolor'] = 'k'
pconf [None,None,'grid_alpha'] = 0.4
pconf [None,None,'grid_linestyle'] = '-'
pconf [None,None,'grid_dashes'] = [1000, 1]
pconf [None,None,'grid_latmax'] = 90
pconf [None,None,'meridians'] = np.arange(-180,180,45)
pconf [None,None,'parallels'] = np.arange(-75,76,15)


###############################################
# Temperatures, potential temperatures
pconf [None,'t','cmap'] = cm.RdYlCyBu2

pconf [None,'pt','cmap'] = cm.RdYlCyBu2

pconf [None,'eqpt','cmap'] = cm.RdYlCyBu2


###############################################
# Moisture
pconf [None,'q','cmap'] = cm.WYlBuKPi
pconf [None,'q','hook'] = lambda q: q*1e3 			# From kg/kg to g/kg
pconfl[None,'q','hook'] = lambda q: q*1e3 			# From kg/kg to g/kg

pconf [None,'tcw','cmap'] = cm.WYlBuKPi
pconf [None,'tcwv','cmap'] = cm.WYlBuKPi


###############################################
# Precipitation
pconf [None,'cp','cmap'] = cm.WBuPi2

pconf [None,'lsp','cmap'] = cm.WBuPi2


###############################################
# Vorticity, divergence
pconf [None,'pv','cmap'] = cm.WYlRd
pconf [None,'pv','hook'] = lambda pv: pv*1e6 		# From SI units to PVU
pconf [None,'pv','scale'] = np.arange(0,11,1)
pconfl[None,'pv','hook'] = lambda pv: pv*1e6 		# From SI units to PVU
pconfl[None,'pv','scale'] = np.array([-2, -1, 1, 2])


###############################################
# Deformation
pconf [None,'defabs','cmap'] = cm.defabs

pconf [None,'defang','scale'] = np.arange(-4,6)*np.pi/8 - np.pi/16 
pconf [None,'defang','colors'] = ['w', (0.625,0.625,0.85), (0.25, 0.25, 0.7), (0.375,0.375,0.6), '0.5', (0.6, 0.6, 0.375), (0.7, 0.7, 0.25), (0.85, 0.86, 0.625), 'w'] 
pconf [None,'defang','ticks'] = np.arange(-4,5)*np.pi/8.0 
pconf [None,'defang','ticklabels'] = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']
pconf [None,'defang','extend'] = 'neither'

pconf [None,'defanr','scale'] = np.arange(-4,6)*np.pi/8 - np.pi/16 
pconf [None,'defanr','colors'] = ['w', (0.625,0.625,0.85), (0.25, 0.25, 0.7), (0.375,0.375,0.6), '0.5', (0.6, 0.6, 0.375), (0.7, 0.7, 0.25), (0.85, 0.86, 0.625), 'w'] 
pconf [None,'defanr','ticks'] = np.arange(-4,5)*np.pi/8.0 
pconf [None,'defanr','ticklabels'] = [u'-π/2', u'-3π/8', u'-π/4', u'-π/8', u'0', u'π/8', u'π/4', u'3π/8', u'π/2']
pconf [None,'defanr','extend'] = 'neither'


###############################################
# Wind speed and direction
pconf [None,'ff','cmap'] = cm.WYlGn
pconf [None,'dd','scale'] = np.arange(0,36.1)*10.0
pconf [None,'dd','ticks'] = np.arange(0,8)*360.0/8.0
pconf [None,'dd','ticklabels'] = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']


###############################################
# Jet and feature detection densities
pconf [None,'jetaxis_freq','hook'] = lambda x: x*1.0e9  # From m/m^2 to km / (1000 km)^2
pconf [None,'jetaxis_freq','scale'] = np.arange(0,3001,250)
pconf [None,'jetaxis_freq','ticks'] = np.arange(0,3001,1000)
pconf [None,'jetaxis_freq','ticklabels'] = list(np.arange(0,3001,1000))
pconf [None,'jetaxis_freq','extend'] = 'max'

pconf [None,'rwb_a','cmap'] = cm.Blues
pconf [None,'rwb_c','cmap'] = cm.Oranges

pconf [None,'block_freq','cmap'] = cm.Blues
pconf [None,'cyc_freq','cmap'] = cm.Blues


###############################################
# Geopotential, sea-level pressure
pconf [None,'z','cmap'] = cm.KPiRdYl1
pconf [None,'z','hook'] = lambda msl: msl/98.1		# From m2s-2 to gpdm
pconfl[None,'z','hook'] = lambda msl: msl/98.1		# From m2s-2 to gpdm

pconf [None,'msl','hook'] = lambda msl: msl/100.0 		# From Pa to hPa
pconfl[None,'msl','hook'] = lambda msl: msl/100.0 		# From Pa to hPa
pconfl[None,'msl','linewidths'] = 0.5
pconfl[None,'msl','scale'] = range(900,1070,5)
pconfl[None,'msl','alpha'] = 0.8


###############################################
# Surface fluxes
shf = [-300,-200,-140,-100,-70,-50,-30,-20,-14,-10, 30,50,70,100,140,200,300,500,700,1000,]
shft = [-10000,]+shf+[10000,]
pconf [None,'sshf','cmap'] = cm.BrBG
pconf [None,'sshf','scale'] = shf
pconf [None,'sshf','norm'] = mpl.colors.BoundaryNorm(shft, 256)
pconf [None,'sshf','extend'] = 'both'
pconf [None,'sshf','cb_tickspacing'] = 'uniform'

pconf [None,'slhf','cmap'] = cm.BrBG
pconf [None,'slhf','scale'] = shf
pconf [None,'slhf','norm'] = mpl.colors.BoundaryNorm(shft, 256)
pconf [None,'slhf','extend'] = 'both'
pconf [None,'slhf','cb_tickspacing'] = 'uniform'


###############################################
# Ororgraphy
pconf [None,'oro','cmap'] = cm.gist_earth
pconf [None,'oro','hook'] = lambda msl: msl/9.81		# From m2s-2 to m
pconf [None,'oro','scale'] = range(-19000,51000,2000)
pconfl[None,'oro','hook'] = lambda msl: msl/9.81		# From m2s-2 to m
pconfl[None,'oro','scale'] = range(10000,80001,10000)


# C'est la fin.
