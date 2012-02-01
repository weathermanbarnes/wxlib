#!/usr/bin/python
#  -*- encoding: utf-8

import datetime
import numpy as np
from metopen import metopen
import static as c
import dynlib

q  = 'defang'
qw = 'defabs'
bench = True

opath = '../deformation/'


for plev in c.plevs:
	nttot = 0
	sum   = np.zeros((151,720))
	sqsum = np.zeros((151,720))
	wsum  = np.zeros((151,720))

	for year in c.years:
		print 'Processing year %d, plev %d' % (year, plev)

		f, dat   = metopen(c.file_std % (year, plev, q), c.q[q])
		if f: f.close()
		f, wgt   = metopen(c.file_std % (year, plev, qw), c.q[qw])
		if f: f.close()

		nt  = dat.shape[0]
		if bench:
			begin = datetime.datetime.now()
		avg, std, wsum = dynlib.stat.basic_weighted(dat, wgt)
		sum  [:,:] += wsum[:,:]*avg[:,:]
		sqsum[:,:] += (nt-1)/nt*wsum[:,:]*std[:,:]**2+(2*nt-1)*avg[:,:]**2
		wsum [:,:] += wsum[:,:]
		nttot      += nt
		if bench:
			print 'Fortran basic stats:', datetime.datetime.now()-begin
		
		ofile = opath+'/'+c.file_stat % (year, plev, q)
		np.savez(ofile, mean=np.ascontiguousarray(avg.astype('f4')), 
				stddev=np.ascontiguousarray(std.astype('f4')) )
		
	
	print 'Saving multi-year stats'
	sum  [:,:]/= nttot
	sqsum[:,:] = np.sqrt((sqsum[:,:]*nttot/wsum[:,:]-(2*nttot-1)*sum[:,:]**2)/(nttot-1))

	ofile = opath+'/'+c.file_mstat % (plev, q)
	np.savez(ofile, mean=np.ascontiguousarray(sum.astype('f4')), 
			stddev=np.ascontiguousarray(sqsum.astype('f4')) )


# the end
