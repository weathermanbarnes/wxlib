#!/usr/bin/python
#  -*- encoding: utf-8

import datetime
import numpy as np
from metopen import metopen
import static as c
import dynlib

q = 'Z'
bench = True

opath = '/scratch/reanalysis'


for plev in c.plevs:
	nttot = 0
	sum   = np.zeros((151,720))
	sqsum = np.zeros((151,720))

	for year in c.years:
		print 'Processing year %d, plev %d' % (year, plev)

		f, dat   = metopen(c.file_std % (year, plev, q), c.q[q])
		nt  = dat.shape[0]
		if bench:
			begin = datetime.datetime.now()
		#avg = dat.mean(axis=0)
		#std = dat.std(axis=0)
		avg, std = dynlib.stat.basic(dat)
		sum  [:,:] += nt*avg[:,:]
		sqsum[:,:] += (nt-1)*std[:,:]**2+(2*nt-1)*avg[:,:]**2
		nttot      += nt
		if bench:
			print 'Fortran basic stats:', datetime.datetime.now()-begin
		
		ofile = opath+'/'+c.file_stat % (year, plev, q)
		np.savez(ofile, mean=avg.astype('f4'), stddev=std.astype('f4'))

		f.close()
	
	print 'Saving multi-year stats'
	sum  [:,:]/= nttot
	sqsum[:,:] = np.sqrt(sqsum[:,:]-(2*nttot-1)*sum[:,:]**2)

	ofile = opath+'/'+c.file_mstat % (plev, q)
	np.savez(ofile, mean=sum.astype('f4'), stddev=sqsum.astype('f4'))


# the end
