#!/usr/bin/env python
#  -*- encoding: utf-8

import os
import numpy as np
import scipy.io.netcdf as nc
import gridlib
import dynlib
import utils
import static as c

from metopen import metopen

#years = range(1979,1988)
#plevs = [1000,950,900,850,800,750,700,650,600,550,500,400,300,200,100]

dt = utils.datetime.datetime

years = c.years
#plevs = c.plevs
plevs = ['pt300', 'pt315', 'pt330', 'pt350']

for year in years:
	for plev in plevs:
		print 'Processing year %d, plev %s' % (year, plev)

		ipath  = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'
		opath  = '/Data/gfi/scratch/csp001/deformation'

		# Open nc file, check if wind data is present
		fu, u  = metopen(c.file_std % (year, plev, 'u'), 'u')
		fv, v  = metopen(c.file_std % (year, plev, 'v'), 'v')

		# Extract grid information 
		grid = gridlib.grid(fu)
		dynlib.config.grid_cyclic_ew = grid.cyclic_ew

		if not u.shape == v.shape:
			raise TypeError, 'Field shape for u wind does not match field shape for v.'

		deff = utils.call(dynlib.diag.def_total, [u,v], grid, cut=c.std_slice, bench=True)

		fu.close()
		fv.close()
		
		ofile = '%s/ei.ans.%d.%s.defabs.npy' % (opath, year, plev)
		begin = dt.now()
		np.save(ofile, np.ascontiguousarray(deff.astype('f4')))
		print 'Saving', dt.now()-begin

		#begin = dt.now()
		#os.spawnl(os.P_WAIT, '/usr/bin/scp', 'scp', ofile, 'gfi063203.klientdrift.uib.no:/media/work/reanalysis/highres')
		#os.spawnl(os.P_WAIT, '/bin/rm', 'rm', ofile)
		#print 'Moving', dt.now()-begin


#
