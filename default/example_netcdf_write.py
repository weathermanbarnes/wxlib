#!/usr/bin/env python
#  -*- encoding: utf-8

# Example script for dynlib:
#  demostrating direct access to reading and writing netcdf files
#
# Example application:
#  shifting the longitude-axis for all chosen netcdf files by 180°,
#  in this example from 0 -- 360° to -180° -- 180°


from dynlib import *
from settings import conf

from scipy.io import netcdf_file as ncf, netcdf_variable as ncv


qs    = ['pv', 'u', 'v']
years = conf.years
plevs = conf.plevs  # for pressure levels
#plevs = conf.ptlevs # for potential temperature levels
#plevs = conf.pvlevs # for potential vorticity levels

# Manual ipath instead of metopen to make sure we get the right nc-file
ipath = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY/'

for year in years:
	for plev in plevs:
		for q in qs:
			print 'Processing', year, plev, q
			
			# Make sure to have set the right ipath and conf.opath!
			inf = ncf(     ipath+conf.file_std % (year, plev, q)+'.nc', 'r')
			otf = ncf(conf.opath+conf.file_std % (year, plev, q)+'.nc', 'w')
			
			# Copying over the basic information of the nc-file
			otf._dims = inf._dims
			otf.dimensions = inf.dimensions
			otf._attributes = inf._attributes
			
			# Copying over the actual variables and possibly adapt them as needed
			for vn, vd in inf.variables.items():
				if vn == q:
					s = vd.data.shape
					dat = np.zeros(s, dtype=vd.data.dtype)
					dat[:,:,:360] = vd.data[:,:,360:]
					dat[:,:,360:] = vd.data[:,:,:360]
					otf.variables[vn] = ncv(dat, vd.typecode(), s, vd.dimensions, vd._attributes)
				elif vn == 'longitude':
					dat = np.arange(-180,180,0.5)
					otf.variables[vn] = ncv(dat, vd.typecode(), vd.data.shape, vd.dimensions, vd._attributes)
				else:
					otf.variables[vn] = vd
			
			# Write the results. DO NOT call both otf.flush() and otf.close() 
			# as this will write all data twice!
			otf.close()

# done
