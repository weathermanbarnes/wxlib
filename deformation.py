#!/usr/bin/python

import numpy as np
import scipy.io.netcdf as nc
import gridlib
import dynlib

ipath  = '/Data/gfi/share/Reanalyses/ERA_INTERIM/DAILY'
opath  = '/scratch/reanalysis'
ncfile = 'ERA_INT_DAILY_wind_700_1979.nc'

# Open nc file, check if wind data is present
f  = nc.netcdf_file('%s/%s' % (ipath, ncfile))
if not 'uwnd' in f.variables or not 'vwnd' in f.variables:
	raise TypeError, 'Expected wind data in netcdf file.'

# Extract grid information 
grid = gridlib.grid(f)
dynlib.config.grid_cyclic_ew = grid.cyclic_ew

# TODO: Check if uwnd and vwnd are on the same grid [see: ['uwnd'].dimensions]
u = f.variables['uwnd'].data
v = f.variables['vwnd'].data

# Calculate deformation
if not grid.nz or grid.nz == 1:
	print '3D mode'
	if len(u.shape) == 4:
		u = u[:,0,:,:]
		v = v[:,0,:,:]
	deff = dynlib.diag.def_angle(u, v, grid.dx, grid.dy)
else:
	print '4D mode'
	raise NotImplementedError
	#for t in len(u.shape[0]):
	#	dylib.diag.def(u[t,:,:,:], v[t,:,:,:], grid.dx, grid.dy)

f.close()

print 'Saving'
np.savez('%s/defang.npz' % opath, defang=deff)

#
