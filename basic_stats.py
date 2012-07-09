#!/usr/bin/python
#  -*- encoding: utf-8

import datetime
import numpy as np
from metopen import metopen
import static as c
import dynlib

q = 'defang'
bench = True

opath = '/work/csp001/deformation'


for plev in c.plevs:
	nttot = 0
	sum   = np.zeros((361,720))
	sqsum = np.zeros((361,720))

	for year in c.years:
		print 'Processing year %d, plev %s' % (year, plev)

		f, dat   = metopen(c.file_std % (year, plev, q), c.q[q])
		nt  = dat.shape[0]
		if bench:
			begin = datetime.datetime.now()
		#avg = dat.mean(axis=0)
		#std = dat.std(axis=0)
		min,max,avg,std = dynlib.stat.basic(dat)
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
		nttot      += nt
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
	sum  [:,:]/= nttot
	sqsum[:,:] = np.sqrt((sqsum[:,:]-(2*nttot-1)*sum[:,:]**2)/(nttot-1))

	ofile = opath+'/'+c.file_mstat % (plev, q)
	np.savez_compressed(ofile, 
		mean   = np.ascontiguousarray(sum.astype('f4')), 
		stddev = np.ascontiguousarray(sqsum.astype('f4')),
		minv   = np.ascontiguousarray(minv.astype('f4')),
		maxv   = np.ascontiguousarray(maxv.astype('f4'))  ) 


# the end
