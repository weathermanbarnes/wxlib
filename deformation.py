#!/usr/bin/python

import numpy as np
import scipy.io.netcdf as nc
import gridlib
import dynlib

ipath  = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'
opath  = '/scratch/reanalysis'
ufile  = 'ei.ans.1979.1000.u.nc'
vfile  = 'ei.ans.1979.1000.v.nc'

# Open nc file, check if wind data is present
fu  = nc.netcdf_file('%s/%s' % (ipath, ufile), 'r')
fv  = nc.netcdf_file('%s/%s' % (ipath, vfile), 'r')
if not 'u' in fu.variables or not 'v' in fv.variables:
	raise TypeError, 'Expected wind data in netcdf file.'

# Extract grid information 
grid = gridlib.grid(fu)
dynlib.config.grid_cyclic_ew = grid.cyclic_ew

u = fu.variables['u']
v = fv.variables['v']
if not u.shape == v.shape:
	raise TypeError, 'Field shape for u wind does not match field shape for v.'

# Calculate deformation
if not grid.nz or grid.nz == 1:
	print '3D mode'
	if len(u.shape) == 4:
		raise NotImplementedError
		u = u[:,0,:,:]
		v = v[:,0,:,:]
	if hasattr(u, 'scale_factor') or hasattr(u, 'add_offset'):
		begin = datetime.datetime.now()
		u_dat = dynlib.conv.scaleoff(u[::], getattr(u, 'scale_factor', 1.0), getattr(u, 'add_offset', 0.0))
		print datetime.datetime.now()-begin
	else:
		u_dat = u[::]
	if hasattr(v, 'scale_factor') or hasattr(v, 'add_offset'):
		#v_dat = dynlib.conv.scaleoff(v[::], getattr(v, 'scale_factor', 1.0), getattr(v, 'add_offset', 0.0))
		begin = datetime.datetime.now()
		v_dat = u[::]*getattr(v, 'scale_factor', 1.0) + getattr(v, 'add_offset', 0.0)
		print datetime.datetime.now()-begin
	else:
		v_dat = v[::]

	deff = dynlib.diag.def_angle(u_dat, v_dat, grid.dx, grid.dy)
else:
	print '4D mode'
	raise NotImplementedError
	#for t in len(u.shape[0]):
	#	dylib.diag.def(u[t,:,:,:], v[t,:,:,:], grid.dx, grid.dy)

fu.close()
fv.close()

print 'Saving'
np.savez('%s/defang.npz' % opath, defang=deff.astype('f4'))

#
