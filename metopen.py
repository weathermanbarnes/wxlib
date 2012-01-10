#!/usr/bin/python

import static
import numpy as np
import os.path

def metopen(filename):
	for path in static.datapath:
		if os.path.exists(path+'/'+filename):
			return np.load(path+'/'+filename)
	
	raise RuntimeError, '%s not found in any data location.' % filename

# the end
