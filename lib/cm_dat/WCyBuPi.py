#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red': ((0.0, 1.0, 1.0),
		   (0.1, 0.3, 0.3),
		   (0.3, 0.1, 0.1),
		   (0.5, 0.1, 0.1),
		   (0.6, 0.1, 0.1),
		   (0.8, 0.8, 0.8),
		   (1.0, 0.4, 0.4)),

	 'green': ((0.0, 1.0, 1.0), 
		   (0.3, 0.9, 0.9),
		   (0.5, 0.65,0.65),
		   (0.6, 0.4, 0.4),
		   (0.8, 0.1, 0.1),
		   (1.0, 0.1, 0.1)),

	 'blue':  ((0.0, 1.0, 1.0),
		   (0.3, 0.7, 0.7),
		   (0.5, 0.8, 0.8),
		   (0.6, 1.0, 1.0),
		   (0.8, 0.7, 0.7),
		   (1.0, 0.4, 0.4))
}

WCyBuPi = LinearSegmentedColormap('WCyBuPi', cm_data)
WCyBuPi_r = LinearSegmentedColormap('WCyBuPi_r', reverse_cm(cm_data))

# Make a clean "from <colormap> import *"
__all__ = ['WCyBuPi', 'WCyBuPi_r']

# For use with viscm
test_cm = WCyBuPi

# the end
