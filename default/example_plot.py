#!/usr/bin/env python
#  -*- encoding: utf-8

# Example script for dynlib:
#  demostrating the use of the plotting functions
#
# Example application:
#  Plot the evolution of winter storm Dagmar:
#  Temperature as shading, wind barbs and geopotential contours on top


from dynlib import *
from settings import *

from datetime import datetime as dt


plotname = 'dagmar'					# Prefix for the all filenames

dates = [dt(2011,12,23,00)+i*td(0.25) for i in range(21)]  # 6-hourly plots for 5 days
plevs = ['500', '700', '850', ] 			# (pressure) levels to plot
proj  = NAmap 						# Use a predefined projection for the North Atlantic


# 
for plev in plevs:
	# Retrieve the data
	T = fig._get_instantaneous('T', dates, plevs=plev, tavg=False)
	u = fig._get_instantaneous('u', dates, plevs=plev, tavg=False)
	v = fig._get_instantaneous('v', dates, plevs=plev, tavg=False)
	
	# On potential temperature levels plot PV contours on top
	if plev[:2] == 'pt':
		con = f._get_instantaneous('pv', dates, plevs=plev, tavg=False)
		Z   = f._get_instantaneous('Z', dates, plevs=plev, tavg=False)
		con_scale = np.array([-2,-1, 1,2])*1e-6

	# On PV levels plot geopotential on top
	elif plev[:2] == 'pv':
		con = f._get_instantaneous('Z', dates, plevs=plev, tavg=False)
		Z   = con
		con_scale = np.arange(1000,100001,1000)
	
	# On presure levels plot geopotential on top
	else:
		con = f._get_instantaneous('Z', dates, plevs=plev, tavg=False)
		Z   = con
		con_scale = np.arange(1000,100001,1000)
	
	# Do the actual plots
	for date, tidx in zip(dates, range(dabs.shape[0])):
		# Figure setup
		fig = plt.figure(figsize=(12.0,8.8), dpi=96)
		fig.subplotpars.update(left=0.01, right=0.99, top=0.99, bottom=0.01)

		# Contours to put on top of the shading and barbs
		overlays = [f.map_overlay_dat(con[tidx,:,:], **conf.contour.merge('default', scale=con_scale) ), ]
		plotconf = conf.contourf.merge('T', overlays=overlays, show=False,
				save='%s/deform_typ_%s_%s_%s.png' % (conf.opath, plotname, plev, date.strftime('%Y%m%d%H'))
		f.map_oro_barb(u[tidx,:,:], v[tidx,:,:], dat=T[tidx,:,:], **plotconf)

	# Clean up
	plt.close('all')


#
