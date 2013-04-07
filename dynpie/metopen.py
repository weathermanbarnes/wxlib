#!/usr/bin/env python

import os.path
import numpy as np
import scipy.io.netcdf as nc
import scipy.io.matlab as mat
from settings import conf as c
import utils 

def metopen(filename, q, cut=c.std_slice, verbose=False, no_dtype_conversion=False):
	for path in c.datapath:
		if verbose:
			print 'Trying: '+path+'/'+filename+'.*'
		if os.path.exists(path+'/'+filename+'.npy'):
			dat = np.load(path+'/'+filename+'.npy', mmap_mode='r')
			dat = dat[cut]
			print 'Found '+path+'/'+filename+'.npy'
			f = None
		elif os.path.exists(path+'/'+filename+'.npz'):
			f   = np.load(path+'/'+filename+'.npz')
			dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.npz'
		elif os.path.exists(path+'/'+filename+'.mat'):
			f   = mat.loadmat(path+'/'+filename+'.mat')
			dat = f[q][cut]
			print 'Found '+path+'/'+filename+'.mat'
		elif os.path.exists(path+'/'+filename+'.nc'):
			f   = nc.netcdf_file(path+'/'+filename+'.nc', 'r')
			var = f.variables[q]
			dat = utils.scale(var, cut)
			print 'Found '+path+'/'+filename+'.nc'
		
		if not no_dtype_conversion:
			dat = dat.astype('f8')

		return f, dat
	
	raise RuntimeError, '%s.npz/.nc not found in any data location.' % filename

# the end
