#!/usr/bin/python
#  -*- encoding: utf-8

import os
import numpy as np
import scipy.io.netcdf as nc
import gridlib
import dynlib
import utils
import static

#years = range(1979,1988)
#plevs = [1000,950,900,850,800,750,700,650,600,550,500,400,300,200,100]

dt = utils.datetime.datetime


for year in static.years:
	for plev in static.plevs:
		print 'Processing year %d, plev %d' % (year, plev)

		ipath  = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'
		opath  = '/scratch/csp001/reanalysis'
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

		u = fu.variables['u']
		v = fv.variables['v']
		if not u.shape == v.shape:
			raise TypeError, 'Field shape for u wind does not match field shape for v.'

		# Cut out    all times    lat > 15°N       all longitudes
		spacetime = (slice(None), slice(210,None), slice(None))
		deff = utils.call(dynlib.diag.def_angle, [u,v], grid, cut=spacetime, bench=True)

		fu.close()
		fv.close()
		
		ofile = '%s/ei.ans.%d.%d.defang.npz' % (opath, year, plev)
		begin = dt.now()
		np.savez(ofile, defang=deff.astype('f4'))
		print 'Saving', dt.now()-begin

		begin = dt.now()
		os.spawnl(os.P_WAIT, '/usr/bin/scp', 'scp', ofile, 'gfi063203.klientdrift.uib.no:/media/work/reanalysis/highres')
		os.spawnl(os.P_WAIT, '/bin/rm', 'rm', ofile)
		print 'Moving', dt.now()-begin


#
