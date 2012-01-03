#!/usr/bin/python

import scipy.io.netcdf as nc

f  = nc.netcdf_file('/Data/gfi/share/Reanalyses/ERA_INTERIM/DAILY/ERA_INT_DAILY_wind_700_1979.nc')

for var in f.variables:
	v = f.variables[var]
	
	# Print name, tuple of dimension names, tuple of dimension lengths
	print var,': ', v.long_name, v.dimensions, v.shape
	
	# Access the data
	dat = v.data
	print 'Min: ', dat.min(), ', mean: ', dat.mean(), ', max: ', dat.max()
	print ''

f.close()

#
