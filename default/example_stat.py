#!/usr/bin/env python
#  -*- encoding: utf-8

# Example script for dynlib:
#  demostrating the use of the statistics functions
#
# Example application:
#  calculate and save some basic statistics (mean, standard deviation, 
#  minimum, maximum and linear trend + the median, the most frequent 
#  value and the complete histogram based on predefined bins) 
#  for all chosen years and vertical levels


from dynlib import *
from settings import conf

from datetime import datetime as dt


bench = True

q = 'defang'
years = conf.years
plevs = ['300', '500', '800']
#plevs = conf.plevs  # for pressure levels
#plevs = conf.ptlevs # for potential temperature levels
#plevs = conf.pvlevs # for potential vorticity levels
s = (361,720) 	# shape of the array

for plev in plevs:
	# initialising multi-year statistics
	nttot = 0 		# total number of time steps
	sum   = np.zeros(s) 	# sum of the values
	sqsum = np.zeros(s) 	# sum of the squares
	tprod = np.zeros(s) 	# sum of timeindex * value

	# initialising the histograms
	if q in conf.bins:
		hist = np.zeros((len(conf.bins[q]),s[0],s[1]))

	for year in years:
		print 'Processing year %d, plev %s' % (year, plev)
		
		# Load the data
		f, dat = metopen(conf.file_std % (year, plev, q), conf.q[q])
		
		# Number of time steps in this file
		nt = dat.shape[0]

		# calculate the basic statistics 
		if bench:
			begin = datetime.datetime.now()
		min,max,avg,std,tprod_out = dynlib.stat.basic(dat, nttot)
		if bench:
			print 'Fortran basic stats:', datetime.datetime.now()-begin
		
		# Save the relevant information for the multi-year statistics
		tprod += tprod_out
		if q in conf.bins:
			mfv,his,med     = dynlib.stat.binned(dat, conf.bins[q])
			hist += his
		sum  [:,:] += nt*avg[:,:]
		sqsum[:,:] += (nt-1)/nt*std[:,:]**2 + avg[:,:]**2
		if year == conf.years[0]:
			minv = min[:,:]
			maxv = max[:,:]
		else:
			minv [:,:]  = np.minimum(min, minv)
			maxv [:,:]  = np.maximum(max, maxv)
		nttot += nt
		
		# Save the yearly statisitcs, make sure to have conf.opath set correctly!
		ofile = conf.opath+'/'+conf.file_stat % (year, plev, q)
		if q in conf.bins:
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
		
		# Close input file
		if f:
			f.close()
	
	# Calculating multi-year statistics
	tmean   =  np.arange(nttot).mean()
	tsqmean = (np.arange(nttot)**2).mean()
	sum  [:,:] /= nttot
	sqsum[:,:] /= len(c.years)
	tprod[:,:] /= nttot
	
	# Trend and y-axis intercept
	trend  = (tprod[:,:] - tmean*sum[:,:]) / (tsqmean - tmean**2)
	icept  = (sum[:,:] - trend[:,:]*tmean)
	
	# Estimated standard deviation of trend and y-axis intercept
	ssqeps = float(nttot)/(nttot-2.0) * (sqsum[:,:] - sum[:,:]**2 - trend[:,:]**2*(tsqmean - tmean**2))
	strend = ssqeps/(nttot*(tsqmean - tmean**2))
	sicept = strend*tsqmean
	
	# Standard deviation from sum of squares
	sqsum[:,:] = np.sqrt((sqsum[:,:] - sum[:,:]**2)*nttot/(nttot-1.0))
	
	# Saving the multi-year statistics
	print 'Saving multi-year stats'
	ofile = opath+'/'+c.file_mstat % (plev, q)
	if q not in c.bins:
		np.savez_compressed(ofile, 
			mean   = np.ascontiguousarray(sum.astype('f4')), 
			stddev = np.ascontiguousarray(sqsum.astype('f4')),
			minv   = np.ascontiguousarray(minv.astype('f4')),
			maxv   = np.ascontiguousarray(maxv.astype('f4')),
			trend  = np.ascontiguousarray(trend.astype('f4')), 
			icept  = np.ascontiguousarray(icept.astype('f4')),
			strend  = np.ascontiguousarray(strend.astype('f4')), 
			sicept  = np.ascontiguousarray(sicept.astype('f4')) )
	else:
		mfv = utils.cal_mfv(hist, conf.bins[q])
		np.savez_compressed(ofile, 
			mean   = np.ascontiguousarray(sum.astype('f4')), 
			stddev = np.ascontiguousarray(sqsum.astype('f4')),
			minv   = np.ascontiguousarray(minv.astype('f4')),
			maxv   = np.ascontiguousarray(maxv.astype('f4')),
			trend  = np.ascontiguousarray(trend.astype('f4')), 
			icept  = np.ascontiguousarray(icept.astype('f4')),
			strend  = np.ascontiguousarray(strend.astype('f4')), 
			sicept  = np.ascontiguousarray(sicept.astype('f4')),
			mfv    = np.ascontiguousarray(mfv.astype('f4')), 
			hist   = np.ascontiguousarray(hist.astype('f4'))  ) 
		


# the end
