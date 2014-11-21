#!/usr/bin/env python
# -*- encoding: utf-8

import matplotlib as mpl

# Make the standard matplotlib colorbars available in the same name space
from matplotlib.cm import *


def greys():
	cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3)),
		 'green': ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3)),
		 'blue':  ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3))  }

	return mpl.colors.LinearSegmentedColormap('my_grey',cdict,256)

def greys_r():
	cdict = {'red':   ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0))  }

	return mpl.colors.LinearSegmentedColormap('my_grey_r',cdict,256)

def greys2():
	cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'green': ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
		 'blue':  ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1))  }

	return mpl.colors.LinearSegmentedColormap('my_grey',cdict,256)


def div_bw():
	cdict = {'red':   ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2)),
		 'green': ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2))  }

	return mpl.colors.LinearSegmentedColormap('my_div_grey',cdict,256)

def div_bw_inv():
	cdict = {'red':   ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0))  }

	return mpl.colors.LinearSegmentedColormap('my_div_grey',cdict,256)


def defabs():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.4, 0.4), (0.867, 1.0, 1.0), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.5, 0.5), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 1.0, 1.0), (0.867, 0.0, 0.0), (1.0, 0.2, 0.2))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def defabs2():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def q():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.30, 0.30),  (0.867, 0.1, 0.1), (1.0, 0.5, 0.5)),
		 'green': ((0.0, 1.0, 1.0), (0.33, 0.65, 0.65),  (0.867, 0.2, 0.2), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 1.0, 1.0), (0.33, 0.80, 0.80),  (0.867, 0.6, 0.6), (1.0, 0.8, 0.8))  }

	return mpl.colors.LinearSegmentedColormap('my_defabs',cdict,256)

def periodic():
	cdict = {'red':   ((0.0, 0.0, 0.0), (0.25, 0.8, 0.8), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0)),
		 'green': ((0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.9, 0.9), (1.0, 0.0, 0.0)),
		 'blue':  ((0.0, 0.6, 0.6), (0.25, 0.0, 0.0), (0.5, 0.2, 0.2), (0.75, 0.0, 0.0), (1.0, 0.6, 0.6))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic',cdict,256)

def periodic2():
	cdict = {'red':   ((0.0, 0.3, 0.3), (0.25, 0.4, 0.4), (0.5, 0.6, 0.6), (0.75, 0.8, 0.8), (1.0, 0.3, 0.3)),
		 'green': ((0.0, 0.2, 0.2), (0.25, 0.8, 0.8), (0.5, 0.6, 0.6), (0.75, 0.4, 0.4), (1.0, 0.2, 0.2)),
		 'blue':  ((0.0, 0.5, 0.5), (0.25, 0.0, 0.0), (0.5, 1.0, 1.0), (0.75, 0.0, 0.0), (1.0, 0.5, 0.5))  }

	return mpl.colors.LinearSegmentedColormap('my_periodic2',cdict,256)


def periodic3():
	cdict = {'red':   ((0.0, 1.0, 1.0), (0.03, 0.9, 0.9), (0.27, 0.3, 0.2), (0.515, 0.5, 0.5), (0.76, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'green': ((0.0, 1.0, 1.0), (0.03, 0.9, 0.9), (0.27, 0.3, 0.2), (0.515, 0.5, 0.5), (0.76, 0.65, 0.7), (1.0, 1.0, 1.0)),
		 'blue':  ((0.0, 1.0, 1.0), (0.03, 1.0, 1.0), (0.27, 0.7, 0.7), (0.515, 0.5, 0.5), (0.76, 0.20, 0.3), (1.0, 0.9, 1.0)) }

	return mpl.colors.LinearSegmentedColormap('my_periodic3',cdict,256)

# that's it
