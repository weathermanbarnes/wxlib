#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red': ((0.0, 1.0, 1.0),
		   (0.33,1.0, 1.0),
		   (0.66, 1.0, 1.0),
		   (1.0, 0.4, 0.4)),

	 'green': ((0.0, 1.0, 1.0),
		   (0.33, 0.8, 0.8),
		   (0.66, 0.2, 0.2), 
		   (1.0, 0.0, 0.0)),

	 'blue':  ((0.0, 1.0, 1.0),
		   (0.33, 0.1, 0.1),
		   (0.66, 0.1, 0.1), 
		   (1.0, 0.2, 0.2)),
}

WYlRd = LinearSegmentedColormap('WYlRd', cm_data)
WYlRd_r = LinearSegmentedColormap('WYlRd_r', reverse_cm(cm_data))

# Make a clean "from <colormap> import *"
__all__ = ['WYlRd', 'WYlRd_r']

# For use with viscm
test_cm = WYlRd

# the end
