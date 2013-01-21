#!/usr/bin/python
#  -*- encoding: utf-8

import os
import numpy as np
import scipy.io.netcdf as nc
import scipy.io as sp
import gridlib
import dynlib
import utils
import static as c


#years = range(1979,1988)
#plevs = [1000,950,900,850,800,750,700,650,600,550,500,400,300,200,100]

dt = utils.datetime.datetime


for year in c.years:
	for plev in c.plevs:
		print 'Processing year %d, plev %s' % (year, plev)

		ipath  = '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY/'
		opath  = './out'
		ufile  = 'ei.ans.%d.%s.u.nc' % (year, plev)
		vfile  = 'ei.ans.%d.%s.v.nc' % (year, plev)
                pvfile = 'ei.ans.%d.%s.pv.nc' % (year, plev)

		# Open nc file, check if wind data is present
		fu  = nc.netcdf_file('%s/%s' % (ipath, ufile), 'r')
		fv  = nc.netcdf_file('%s/%s' % (ipath, vfile), 'r')
                fpv = nc.netcdf_file('%s/%s' % (ipath, pvfile), 'r')
		if not 'u' in fu.variables or not 'v' in fv.variables or not 'pv' in fpv.variables:
			raise TypeError, 'Expected data not in netcdf file.'

		# Extract grid information 
		grid = gridlib.grid(fu)
		dynlib.config.grid_cyclic_ew = grid.cyclic_ew

		u = fu.variables['u']
		v = fv.variables['v']
                pv = fpv.variables['pv']
		if not u.shape == v.shape:
			raise TypeError, 'Field shape for u wind does not match field shape for v.'
                if not u.shape == pv.shape:
                        raise TypeError, 'Field shape for pv does not match field shape for u.'
		#deff = utils.call(dynlib.diag.def_angle, [u,v], grid, cut=c.std_slice, bench=True)

		fu.close()
		fv.close()
                fpv.close()
		
                deftot = utils.call(dynlib.diag.def_total, [u,v], grid, cut=c.std_slice, bench=True)
                [pvgradx,pvgrady] = utils.call(dynlib.diag.grad, [pv], grid, cut=c.std_slice, bench=True)
                pvgrad=np.sqrt(pvgradx*pvgradx+pvgrady*pvgrady)
                divergence=utils.call(dynlib.diag.div, [u,v], grid, cut=c.std_slice, bench=True)
                ofile = '%s/ei.ans.%d.%s.shear.mat' % (opath, year, plev)
                #print u.typecode() # is h
                udat=utils.scale(u, cut=c.std_slice, bench=False);
                vdat=utils.scale(v, cut=c.std_slice, bench=False); 
                pvdat=utils.scale(pv, cut=c.std_slice, bench=False);
                matdict = {'u':udat, 'v':vdat, 'pv':pvdat, 'deftot': deftot, 'pvgrad':pvgrad,'divergence':divergence }
                #udat2=udat.view(dtype=np.float64)
                sp.savemat(ofile,matdict, format='5')
                #sp.savemat(ofile,{'u':udat},format='5')

                # save .npz
		#ofile = '%s/ei.ans.%d.%s.defang.npz' % (opath, year, plev)
		#begin = dt.now()
		#np.savez(ofile, defang=deff.astype('f4'))
		#print 'Saving', dt.now()-begin

		#begin = dt.now()
		#os.spawnl(os.P_WAIT, '/usr/bin/scp', 'scp', ofile, 'gfi063203.klientdrift.uib.no:/media/work/reanalysis/highres')
		#os.spawnl(os.P_WAIT, '/bin/rm', 'rm', ofile)
		#print 'Moving', dt.now()-begin


#
