#!/usr/bin/python
#  -*- encoding: utf-8

import numpy as np
from scipy.io import netcdf_file as ncf, netcdf_variable as ncv
import static as c

years = [2011, ] #c.years
plevs = ['pt320', ] #['pt300', 'pt315', 'pt330', 'pt350']
qs    = ['pv', ] #['pv', 'u', 'v']

ipath = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY/'
opath = '/Data/gfi/work/csp001/pt_lonshift/'

for year in years:
	for plev in plevs:
		for q in qs:
			print year, plev, q
			inf = ncf(ipath+c.file_std % (year, plev, q)+'.nc', 'r')
			otf = ncf(opath+c.file_std % (year, plev, q)+'.nc', 'w')

			otf._dims = inf._dims
			otf.dimensions = inf.dimensions
			otf._attributes = inf._attributes

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

			otf.close()

# done
