#!/usr/bin/python
#  -*- encoding: utf-8

import os
import numpy as np
import scipy.io.netcdf as nc
import gridlib
import dynlib
import utils
import static as c

#years = range(1979,1988)
#plevs = [1000,950,900,850,800,750,700,650,600,550,500,400,300,200,100]

dt = utils.datetime.datetime


for year in [1979, ]: #c.years:
	for plev in c.plevs:
		print 'Processing year %d, plev %d' % (year, plev)

		ipath  = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'
		opath  = '/Data/gfi/scratch/csp001/deformation'
		ufile  = 'ei.ans.%d.%d.u.nc' % (year, plev)
		vfile  = 'ei.ans.%d.%d.v.nc' % (year, plev)

		# Open nc file, check if wind data is present
		fu  = nc.netcdf_file('%s/%s' % (ipath, ufile), 'r')
		fv  = nc.netcdf_file('%s/%s' % (ipath, vfile), 'r')
		if not 'u' in fu.variables or not 'v' in fv.variables:
			raise TypeError, 'Expected wind data in netcdf file.'

		# Extract grid information 
		grid = gridlib.grid(fu)
		dynlib.config.grid_cyclic_ew = grid.cyclic_ew

		u = fu.variables['u'].data[::]
		v = fv.variables['v'].data[::]
		if not u.shape == v.shape:
			raise TypeError, 'Field shape for u wind does not match field shape for v.'

		deff = utils.call(dynlib.diag.def_angle, [u,v], grid, cut=c.std_slice, bench=True)

		fu.close()
		fv.close()
		
		ofile = '%s/ei.ans.%d.%d.defang.npy' % (opath, year, plev)
		begin = dt.now()
		np.save(ofile, np.ascontiguousarray(deff.astype('f4')))
		print 'Saving', dt.now()-begin

		#begin = dt.now()
		#os.spawnl(os.P_WAIT, '/usr/bin/scp', 'scp', ofile, 'gfi063203.klientdrift.uib.no:/media/work/reanalysis/highres')
		#os.spawnl(os.P_WAIT, '/bin/rm', 'rm', ofile)
		#print 'Moving', dt.now()-begin


#
