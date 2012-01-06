#!/usr/bin/python

import math
import datetime

#
# Automatic scaling according to netcdf attributes "scale_factor" and "add_offset"
def scale(var, cut, bench=False):
	if hasattr(var, 'scale_factor') or hasattr(var, 'add_offset'):
		if bench:
			begin = datetime.datetime.now()
		# Python/numpy version is faster than Fortran function
		var_dat = var[cut]*getattr(var, 'scale_factor', 1.0) + getattr(var, 'add_offset', 0.0)
		#u_dat = dynlib.conv.scaleoff(u_dat, getattr(u, 'scale_factor', 1.0), getattr(u, 'add_offset', 0.0))
		if bench:
			print 'Python scaleoff', datetime.datetime.now()-begin
	
	return var_dat

#
# Generic calculation preparations and the actual call of the Fortran function
def call(func, vars, grid, cut=(slice(None),slice(None),slice(None)), bench=False):
	if not grid.nz or grid.nz == 1:
		print '3D mode'
		if len(vars[0].shape) == 4:
			raise NotImplementedError
		
		args = []
		for var in vars:
			args.append(scale(var, cut, bench=bench))
		
		args.extend([grid.dx[cut[1:]], grid.dy[cut[1:]]])
		if bench:
			begin = datetime.datetime.now()
		deff = func(*args) 
		if bench:
			print 'Calculation', datetime.datetime.now()-begin
	
	elif not grid.nt or grid.nt == 1:
		print '3D mode'
		raise NotImplementedError
	
	else:
		print '4D mode'
		raise NotImplementedError
		#for t in len(u.shape[0]):
		#	dylib.diag.def(u[t,:,:,:], v[t,:,:,:], grid.dx, grid.dy)
	
	return deff

#
# Reimplementation of the recpective function in dynlib.diag for benchmarking.
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


#
