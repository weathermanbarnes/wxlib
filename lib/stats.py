#!/usr/bin/env python

import os
import sys
import pickle

import numpy as np
import dynfor


# Take over the contents of a fortran module
pth = os.path.abspath(os.path.join(os.path.dirname(__file__), '.dynfor_doc.pickle'))
if os.path.isfile(pth):
	f = open(pth, 'r')
	fortran_doc = pickle.load(f)
	f.close()
else:
	print 'Warning: Fortran documentation not found. Make sure you have compiled dynlib.'
	fortran_doc = {}

for name, value in dynfor.stat.__dict__.items():
	# Backup the f2py generated docstring
	value.__f2pydoc__ = value.__doc__
	# Set the doctring from the fortran_doc extracted documentation
	value.__doc__ = fortran_doc.get('stats.%s' % name, None)
	# Make available as a member of this python module
	locals()[name] = value




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
