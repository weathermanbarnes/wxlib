#!/usr/bin/python
#  -*- encoding: utf-8

import os
import numpy as np
import scipy.io.netcdf as nc
import gridlib
import dynlib
import utils
import static as c

from metopen import metopen

dt = utils.datetime.datetime

years = [1979, ] #c.years
#plevs = c.plevs
plevs = ['pt300', 'pt315', 'pt330', 'pt350']

for year in years:
	for plev in plevs:
		print 'Processing year %d, plev %s' % (year, plev)

		opath  = '/Data/gfi/work/csp001/geop_from_montgp'

		# Open nc file, check if wind data is present
		begin = dt.now()
		fm, m  = metopen(c.file_std % (year, plev, 'm'), c.q['m'])
		fp, p  = metopen(c.file_std % (year, plev, 'p'), c.q['p'])
		print 'Loading', dt.now()-begin

		# Extract grid information 
		grid = gridlib.grid(fm)
		dynlib.config.grid_cyclic_ew = grid.cyclic_ew

		if not m.shape == p.shape:
			raise TypeError, 'Field shape for m does not match field shape for p.'

		theta = np.ones(m.shape) * int(plev[2:])

		res = utils.call(dynlib.diag.geop_from_montgp, [m,theta,p], grid, cut=c.std_slice, bench=True)
		res, scale, off = utils.unscale(res)

		#ofile = '%s/ei.ans.%d.%s.defabs.npy' % (opath, year, plev)
		#begin = dt.now()
		#np.save(ofile, np.ascontiguousarray(res.astype('f4')))
		#print 'Saving', dt.now()-begin

		from scipy.io import netcdf_file as ncf, netcdf_variable as ncv

		ofile = '%s/ei.ans.%d.%s.Z.nc' % (opath, year, plev)
		begin = dt.now()
		fo = ncf(ofile, 'w')
		fo._dims = fm._dims
		fo.dimensions = fm.dimensions
		fo._attributes = fm._attributes
		for vn, vd in fm.variables.items():
			if not vn == c.q['m']:
				fo.variables[vn] = vd
		fo.variables['z'] = ncv(res, 'h', res.shape, fm.variables[c.q['m']].dimensions,
			{'units': 'm**2 s**-2', 'long_name': 'Geopotential', 'add_offset': off, 'scale_factor': scale})
		fo.close()
		print 'Saving', dt.now()-begin

		fm.close()
		fp.close()
		


#
