#!/usr/bin/python

import os.path
import numpy as np
import scipy.io.netcdf as nc
import static as c
import utils 

def metopen(filename, q):
	for path in c.datapath:
		if os.path.exists(path+'/'+filename+'.npz'):
			f   = np.load(path+'/'+filename+'.npz')
			dat = f[q].astype('f8')
			return f, dat
		elif os.path.exists(path+'/'+filename+'.nc'):
			f   = nc.netcdf_file(path+'/'+filename+'.nc', 'r')
			var = f.variables[q]
			dat = utils.scale(var, cut=c.std_slice)
			return f, dat
	
	raise RuntimeError, '%s.npz/.nc not found in any data location.' % filename

# the end
