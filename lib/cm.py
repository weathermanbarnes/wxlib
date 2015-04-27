#!/usr/bin/env python
# -*- encoding: utf-8

import matplotlib as mpl

# Make the standard matplotlib colorbars available in the same name space
from matplotlib.cm import *


def greys():
	''' Colorbar from white to 30% grey 

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	
	See also
	--------
	:meth:`greys_r`
	'''
	cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3)),
		 'green': ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3)),
		 'blue':  ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3))  }

	return mpl.colors.LinearSegmentedColormap('my_grey',cdict,256)

def greys_r():
	''' Colorbar from 30% grey to white

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	
	See also
	--------
	:meth:`greys`
	'''
	cdict = {'red':   ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0))  }

	return mpl.colors.LinearSegmentedColormap('my_grey_r',cdict,256)

def greys2():
	''' Colorbar from white to 10% grey 

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	'''	
	cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'green': ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'blue':  ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1))  }

	return mpl.colors.LinearSegmentedColormap('my_darkgrey',cdict,256)


def div_bw():
	''' Colorbar from 60% grey to white to 20% grey

	Divergent colorbar for black-white figures.

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	
	See also
	--------
	:meth:`div_bw_inv`
	'''	
	cdict = {'red':   ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2)),
		 'green': ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2))  }

	return mpl.colors.LinearSegmentedColormap('my_div_grey',cdict,256)

def div_bw_r():
	''' Colorbar from 60% grey to 20% grey to white

	Divergent colorbar for black-white figures.

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	
	See also
	--------
	:meth:`div_bw`
	'''	
	cdict = {'red':   ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0))  }

	return mpl.colors.LinearSegmentedColormap('my_div_grey_r',cdict,256)


def defabs():
	''' Colorbar from white to 15% grey to red

	Often used for total deformation and jet axis frequencies

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	'''
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)


def q():
	''' Colorbar from white to dark blue to violet

	Often used for humidity, rain measures, etc.

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	'''
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.30, 0.30),  (0.867, 0.1, 0.1), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.65, 0.65),  (0.867, 0.2, 0.2), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 0.80, 0.80),  (0.867, 0.6, 0.6), (1.0, 0.8, 0.8))  }

	return mpl.colors.LinearSegmentedColormap('my_q',cdict,256)

def periodic():
	''' Colorbar from white to blue to 50% grey to yellow back to white

	Useful for periodic variablies like the deformation angle or the wind direction.

	Returns
	-------
	mpl.colors.LinearSegmentedColormap
	    Colormap
	'''
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.25, 0.3, 0.2), (0.5, 0.5, 0.5), (0.75, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.25, 0.3, 0.2), (0.5, 0.5, 0.5), (0.75, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.25, 0.7, 0.7), (0.5, 0.5, 0.5), (0.75, 0.20, 0.3), (1.0, 1.0, 1.0)) }

	return mpl.colors.LinearSegmentedColormap('my_periodic',cdict,256)

# that's it
