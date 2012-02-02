#!/usr/bin/python

import os.path
import numpy as np
import scipy.io.netcdf as nc
import scipy.io.matlab as mat
import static as c
import utils 

def metopen(filename, q, cut=c.std_slice):
	for path in c.datapath:
		if os.path.exists(path+'/'+filename+'.npy'):
			dat = np.load(path+'/'+filename+'.npy', mmap_mode='r')
			dat = dat[cut].astype('f8')
			print 'Found '+path+'/'+filename+'.npy'
			return None, dat
		elif os.path.exists(path+'/'+filename+'.npz'):
			f   = np.load(path+'/'+filename+'.npz')
			dat = f[q].astype('f8')
			print 'Found '+path+'/'+filename+'.npz'
			return f, dat
		elif os.path.exists(path+'/'+filename+'.mat'):
			f   = mat.loadmat(path+'/'+filename+'.mat')
			dat = f[q].astype('f8')
			print 'Found '+path+'/'+filename+'.mat'
			return f, dat
		elif os.path.exists(path+'/'+filename+'.nc'):
			f   = nc.netcdf_file(path+'/'+filename+'.nc', 'r')
			var = f.variables[q]
			dat = utils.scale(var, cut=c.std_slice)
			print 'Found '+path+'/'+filename+'.nc'
			return f, dat
	
	raise RuntimeError, '%s.npz/.nc not found in any data location.' % filename

# the end
