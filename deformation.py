#!/usr/bin/python

import numpy as np
import scipy.io.netcdf as nc
import gridlib
import dynlib
import datetime
import math

ipath  = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'
opath  = '/scratch/csp001/reanalysis'
ufile  = 'ei.ans.1979.700.u.nc'
vfile  = 'ei.ans.1979.700.v.nc'

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
begin = datetime.datetime.now()
u_dat = u[::]
v_dat = v[::]
print 'Reading in', datetime.datetime.now()-begin

def def_angle(u_dat, v_dat, grid):
	deff = np.zeros(u_dat.shape)
	for k in range(u_dat.shape[0]):
		for j in range(1,grid.ny-1):
			for i in range(1,grid.nx-1):
				def_shear   = (u[k,j+1,i]-u[k,j-1,i])/grid.dy[j,i] \
					    + (v[k,j,i+1]-v[k,j,i-1])/grid.dx[j,i]
				def_stretch = (u[k,j,i+1]-u[k,j,i-1])/grid.dx[j,i] \
					    - (v[k,j+1,i]-v[k,j-1,i])/grid.dy[j,i]
				deff[k,j,i] = 0.5*math.atan2(def_shear, def_stretch)
			if grid.cyclic_ew:
				def_shear   = (u[k,j+1,i]-u[k,j-1,i])/grid.dy[j,1] \
					    + (v[k,j,  1]-u[k,j, -1])/grid.dx[j,1]
				def_stretch = (u[k,j,  1]-u[k,j, -1])/grid.dx[j,1] \
					    - (v[k,j+1,i]-v[k,j-1,i])/grid.dy[j,1]
				deff[k,j,1] = 0.5*math.atan2(def_shear, def_stretch)
				def_shear   = (u[k,j+1,i]-u[k,j-1,i])/grid.dy[j,1] \
					    + (v[k,j,  0]-u[k,j, -2])/grid.dx[j,1]
				def_stretch = (u[k,j,  0]-u[k,j, -2])/grid.dx[j,1] \
					    - (v[k,j+1,i]-v[k,j-1,i])/grid.dy[j,1]
				deff[k,j,-1] = 0.5*math.atan2(def_shear, def_stretch)

	return deff

# Calculate deformation
if not grid.nz or grid.nz == 1:
	print '3D mode'
	if len(u.shape) == 4:
		raise NotImplementedError
		u = u[:,0,:,:]
		v = v[:,0,:,:]
	if hasattr(u, 'scale_factor') or hasattr(u, 'add_offset'):
		begin = datetime.datetime.now()
		#u_dat = dynlib.conv.scaleoff(u_dat, getattr(u, 'scale_factor', 1.0), getattr(u, 'add_offset', 0.0))
		u_dat = u_dat*getattr(u, 'scale_factor', 1.0) + getattr(u, 'add_offset', 0.0)
		print 'Python scaleoff', datetime.datetime.now()-begin
	
	if hasattr(v, 'scale_factor') or hasattr(v, 'add_offset'):
		begin = datetime.datetime.now()
		v_dat = dynlib.conv.scaleoff(v_dat, getattr(v, 'scale_factor', 1.0), getattr(v, 'add_offset', 0.0))
		#v_dat = v_dat*getattr(v, 'scale_factor', 1.0) + getattr(v, 'add_offset', 0.0)
		print 'Fortran scaleoff', datetime.datetime.now()-begin
	
	begin = datetime.datetime.now()
	deff = dynlib.diag.def_angle(u_dat, v_dat, grid.dx, grid.dy)
	print 'Fortran deformation', datetime.datetime.now()-begin
	#begin = datetime.datetime.now()
	#deff = def_angle(u_dat, v_dat, grid)
	#print 'Python deformation', datetime.datetime.now()-begin

	
else:
	print '4D mode'
	raise NotImplementedError
	#for t in len(u.shape[0]):
	#	dylib.diag.def(u[t,:,:,:], v[t,:,:,:], grid.dx, grid.dy)

fu.close()
fv.close()

begin = datetime.datetime.now()
np.savez('%s/defang.1979.700.npz' % opath, defang=deff.astype('f4'))
print 'Saving', datetime.datetime.now()-begin

#
