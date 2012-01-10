#!/usr/bin/python
#  -*- encoding: utf-8

import datetime
import numpy as np
from metopen import metopen
import static
import dynlib

q = 'defang'
bench = True

opath = '/scratch/reanalysis'


for plev in static.plevs:
	nttot = 0
	sum   = np.zeros((151,720))
	sqsum = np.zeros((151,720))

	for year in static.years:
		print 'Processing year %d, plev %d' % (year, plev)

		f   = metopen('ei.ans.%d.%d.%s.npz' % (year, plev, q))
		dat = f[q].astype('f8')
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
		
		ofile = '%s/ei.ans.%d.%d.%s.stats.npz' % (opath, year, plev, q)
		np.savez(ofile, mean=avg.astype('f4'), stddev=std.astype('f4'))
	
	print 'Saving multi-year stats'
	sum  [:,:]/= nttot
	sqsum[:,:] = np.sqrt(sqsum[:,:]-(2*nttot-1)*sum[:,:]**2)

	ofile = '%s/ei.ans.mean.%d.%s.stats.npz' % (opath, plev, q)
	np.savez(ofile, mean=sum.astype('f4'), stddev=sqsum.astype('f4'))


# the end
