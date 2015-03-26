#!/usr/bin/env python

import sys

import dynfor
import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.stat, 'stat', sys.modules[__name__])


import numpy as np



def running_mean(ts, rlen=8):
	rmean = np.zeros((len(ts)-rlen+1,))
	for i in range(len(ts)-rlen+1):
		rmean[i] = sum(ts[i:i+rlen])/rlen
	
	return rmean


def running_mean_periodic(ts, rlen=3):
	rmean = np.zeros(len(ts))
	for i in range(len(ts)):
		for j in np.arange(rlen)-rlen/2:
			rmean[i] += ts[(i+j)%len(ts)]
	
	rmean /= rlen

	return rmean


def filter7p(ts):
	rmean = np.zeros((len(ts)-6))
	for i in range(len(ts)-6):
		rmean[i] = (ts[i]-6*ts[i+1]+15*ts[i+2]+44*ts[i+3]+15*ts[i+4]-6*ts[i+5]+ts[i+6])/64.0
	
	return rmean


#
