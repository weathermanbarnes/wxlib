#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red': ((0.0, 0.0, 0.0),
		   (0.15, 0.1, 0.1),
		   (0.35, 0.1, 0.1),
		   (0.5, 1.0, 1.0),
		   (0.65,1.0, 1.0),
		   (0.85, 0.9, 0.9),
		   (1.0, 0.4, 1.0)),

	 'green': ((0.0, 0.0, 0.0),
		   (0.15, 0.2, 0.2),
		   (0.35, 0.8, 0.8),
		   (0.5, 1.0, 1.0),
		   (0.65, 0.8, 0.8),
		   (0.85,0.2, 0.2),
		   (1.0, 0.0, 0.0)),

	 'blue':  ((0.0, 0.2, 0.2),
		   (0.15, 0.6, 0.6),
		   (0.35, 0.8, 0.8),
		   (0.5, 1.0, 1.0),
		   (0.65, 0.1, 0.1),
		   (0.85, 0.1, 0.1),
		   (1.0, 0.0, 0.0))
}

RdYlCyBu2 = LinearSegmentedColormap('RdYlCyBu2', cm_data)
RdYlCyBu2_r = LinearSegmentedColormap('RdYlCyBu2_r', reverse_cm(cm_data))

# Make a clean "from <colormap> import *"
__all__ = ['RdYlCyBu2', 'RdYlCyBu2_r']

# For use with viscm
test_cm = RdYlCyBu2

# the end
