#!/usr/bin/python
#  -*- encoding: utf-8

import datetime
import numpy as np
from metopen import metopen
import static as c
import dynlib

q = 'defabs'
bench = True

opath = '/work/csp001/deformation'


for plev in c.plevs:
	nttot = 0
	sum   = np.zeros((361,720))
	sqsum = np.zeros((361,720))
	tprod = np.zeros((361,720))

	for year in c.years:
		print 'Processing year %d, plev %s' % (year, plev)

		f, dat   = metopen(c.file_std % (year, plev, q), c.q[q])
		nt  = dat.shape[0]
		if bench:
			begin = datetime.datetime.now()
		min,max,avg,std,tprod_out = dynlib.stat.basic(dat, nttot)
		tprod += tprod_out
		if q in c.bins:
			mfv,his,med     = dynlib.stat.binned(dat, c.bins[q])
		sum  [:,:] += nt*avg[:,:]
		sqsum[:,:] += (nt-1)*std[:,:]**2+(2*nt-1)*avg[:,:]**2
		if year == c.years[0]:
			minv = min[:,:]
			maxv = max[:,:]
		else:
			minv [:,:]  = np.minimum(min, minv)
			maxv [:,:]  = np.maximum(max, maxv)
		nttot += nt
		if bench:
			print 'Fortran basic stats:', datetime.datetime.now()-begin
		
		ofile = opath+'/'+c.file_stat % (year, plev, q)
		if q in c.bins:
			np.savez_compressed(ofile, 
				mean   = np.ascontiguousarray(avg.astype('f4')), 
				stddev = np.ascontiguousarray(std.astype('f4')),
				minv   = np.ascontiguousarray(min.astype('f4')),
				maxv   = np.ascontiguousarray(max.astype('f4')), 
				mfv    = np.ascontiguousarray(mfv.astype('f4')), 
				hist   = np.ascontiguousarray(his.astype('f4')), 
				median = np.ascontiguousarray(med.astype('f4'))  )
		else:
			np.savez_compressed(ofile, 
				mean   = np.ascontiguousarray(avg.astype('f4')), 
				stddev = np.ascontiguousarray(std.astype('f4')),
				minv   = np.ascontiguousarray(min.astype('f4')),
				maxv   = np.ascontiguousarray(max.astype('f4'))  )
		
		if f:
			f.close()
	
	print 'Saving multi-year stats'
	tsum   =  np.arange(nttot).sum()
	tsqsum = (np.arange(nttot)**2).sum()

	trend = (tprod[:,:] - 1.0/nttot*(tsum*sum[:,:])) / (tsqsum - 1.0/nttot*tsum**2)
	icept = (sum[:,:] - trend[:,:]*tsum)/nttot

	sum  [:,:]/= nttot
	sqsum[:,:] = np.sqrt((sqsum[:,:]-(2*nttot-1)*sum[:,:]**2)/(nttot-1))

	ofile = opath+'/'+c.file_mstat % (plev, q)
	np.savez_compressed(ofile, 
		mean   = np.ascontiguousarray(sum.astype('f4')), 
		stddev = np.ascontiguousarray(sqsum.astype('f4')),
		minv   = np.ascontiguousarray(minv.astype('f4')),
		maxv   = np.ascontiguousarray(maxv.astype('f4')),
		trend  = np.ascontiguousarray(trend.astype('f4')), 
		icept  = np.ascontiguousarray(icept.astype('f4'))  ) 


# the end
