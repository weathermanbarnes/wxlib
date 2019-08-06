#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red': ((0.0, 1.0, 1.0),
		   (0.4, 0.3, 0.3),
		   (0.9, 0.1, 0.1),
		   (1.0, 0.8, 0.8)),

	 'green': ((0.0, 1.0, 1.0), 
		   (0.4, 0.65, 0.65),
		   (0.9, 0.2, 0.2),
		   (1.0, 0.1, 0.1)),

	 'blue':  ((0.0, 1.0, 1.0),
		   (0.4, 0.8, 0.8),
		   (0.9, 0.4, 0.4),
		   (1.0, 0.7, 0.7))
}

WBuPi2 = LinearSegmentedColormap('WBuPi2', cm_data)
WBuPi2_r = LinearSegmentedColormap('WBuPi2_r', reverse_cm(cm_data))

# Make a clean "from <colormap> import *"
__all__ = ['WBuPi2', 'WBuPi2_r']

# For use with viscm
test_cm = WBuPi2

# the end
