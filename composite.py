#!/usr/bin/python
# -*- encoding: utf-8

import sys
import numpy as np
import figures as f
import static as c
from metopen import metopen

# --------------------------------------------------------------

# Data base
year = 1982
plev = 800

# Deformation events for which point?
yidx = 57
xidx = 371

# Critera for composite
dabs_thres = 0.5/3600.0

# What to do?
PLOT_THRESCNT = False
PLOT_Z = True
PLOT_BARB = False
PLOT_QUIVER = False
PLOT_DEFANG = True

# --------------------------------------------------------------

fabs, dabs = metopen(c.file_std % (year, plev, 'defabs'), c.q['defabs'])
if fabs: fabs.close()

if PLOT_THRESCNT:
	t = dabs[:,:,:] > dabs_thres
	t = t.sum(axis=0)
	t[:3,:]  = np.ma.masked
	t[-4:,:] = np.ma.masked
	f.wmap_oro_dat(t)
	sys.exit(0)

if PLOT_Z or PLOT_BARB or PLOT_QUIVER or PLOT_DEFANG:
	t = dabs[:,yidx,xidx] > dabs_thres
	n = t.sum()
	print 'Looking at %.2f°N, %.2f°E' % (f.lat[yidx,0], f.lon[0,xidx])
	print 'Found %d events, corresponding to %.2f%%' % (n, n/14.6)
	if n < 1:
		print 'Stoping here. Nothing to show.'
		sys.exit(1)

if PLOT_Z:
	fZ, dZ = metopen(c.file_std % (year, plev, 'Z'), c.q['Z'])
	Zmean = dZ.mean(axis=0)
	if fZ : fZ.close()
	dZ = dZ[t,:,:].mean(axis=0)
	f.map_oro_dat(dZ-Zmean, mark=(yidx,xidx), plev=plev)

if PLOT_BARB or PLOT_QUIVER:
	fu, du = metopen(c.file_std % (year, plev, 'u'), c.q['u'])
	fv, dv = metopen(c.file_std % (year, plev, 'v'), c.q['v'])
	umean = du.mean(axis=0)
	vmean = dv.mean(axis=0)
	if fu : fu.close()
	if fv : fv.close()
	du = du[t,:,:].mean(axis=0)
	dv = dv[t,:,:].mean(axis=0)
	dabsmean = dabs[t,:,:].mean(axis=0)
	f.map_oro_barb(du-umean, dv-vmean, dat=dabsmean, mark=(yidx,xidx), quiver=PLOT_QUIVER, plev=plev)

if PLOT_DEFANG:
	fang, dang = metopen(c.file_std % (year, plev, 'defang'), c.q['defang'])
	if fang : fang.close()
	dangmean = dang[t,:,:].mean(axis=0)
	dabsmean = dabs[t,:,:].mean(axis=0)
	f.map_oro_deform(dabsmean, dangmean, mark=(yidx,xidx), plev=plev)



# the end
