#!/usr/bin/env python
#  -*- encoding: utf-8

# Example script for dynlib:
#  demostrating the use of the diagnostics functions
#
# Example application:
#  calculate and save total deformation for all chosen years 
#  and vertical levels


from dynlib import *
from settings import conf

from datetime import datetime as dt


years = conf.years
plevs = conf.plevs  # for pressure levels
#plevs = conf.ptlevs # for potential temperature levels
#plevs = conf.pvlevs # for potential vorticity levels

for year in years:
	for plev in plevs:
		print 'Processing year %d, plev %s' % (year, plev)

		# Open nc file, check if wind data is present
		fu, u, grid  = metopen(conf.file_std % (year, plev, 'u'), 'u')
		fv, v, grid  = metopen(conf.file_std % (year, plev, 'v'), 'v')

		# Extract grid information 
		dynlib.config.grid_cyclic_ew = grid.cyclic_ew
		
		# Consistency check
		if not u.shape == v.shape:
			raise TypeError, 'Field shape for u wind does not match field shape for v.'
		
		# Calculate the actual diagnostc, here def_total
		deff = utils.call(dynlib.diag.def_total, [u,v], grid, cut=conf.std_slice, bench=True)
		
		# Close the file files
		fu.close()
		fv.close()

		# Write results
		# OBS: Remember to set conf.opath to a sensible location in your project settings.py. 
		#      Default is current directory!
		ofile = '%s/ei.ans.%d.%s.defabs.npy' % (conf.opath, year, plev)
		begin = dt.now()
		np.save(ofile, np.ascontiguousarray(deff.astype('f4')))
		print 'Saving', dt.now()-begin


#
