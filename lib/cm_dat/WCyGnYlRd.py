#!/usr/bin/env python
# -*- encoding: utf-8

from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red': ((0.0, 1.0, 1.0),
		   (0.125, 0.4, 0.4),
		   (0.25, 0.1, 0.1),
		   (0.4, 0.9, 0.9),
		   (0.8, 1.0, 1.0),
		   (1.0, 0.5, 0.5)),

	 'green': ((0.0, 1.0, 1.0),
		   (0.1, 1.0, 1.0),
		   (0.25, 0.7, 0.7), 
		   (0.4, 0.9, 0.9),
		   (0.8, 0.1, 0.1),
		   (1.0, 0.1, 0.1)),

	 'blue':  ((0.0, 1.0, 1.0),
		   (0.1, 1.0, 1.0),
		   (0.25, 0.1, 0.1),
		   (0.4, 0.1, 0.1),
		   (0.8, 0.1, 0.1),
		   (1.0, 0.2, 0.2))
}

WCyGnYlRd = LinearSegmentedColormap('WCyGnYlRd', cm_data)
WCyGnYlRd_r = LinearSegmentedColormap('WCyGnYlRd_r', cm_data)

# Make a clean "from <colormap> import *"
__all__ = ['WCyGnYlRd', 'WCyGnYlRd_r']

# For use with viscm
test_cm = WCyGnYlRd

# the end
