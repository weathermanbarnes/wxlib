#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red':   ((0.0, 1.0, 1.0), (0.25, 1.0, 1.0),  (0.65, 0.2, 0.2), (1.0, 0.0, 0.0)),
	   'green': ((0.0, 1.0, 1.0), (0.25, 1.0, 1.0),  (0.65, 0.8, 0.8), (1.0, 0.3, 0.3)),
	   'blue':  ((0.0, 1.0, 1.0), (0.25, 0.2, 0.2),  (0.65, 0.2, 0.2), (1.0, 0.0, 0.0))  }

WYlGn = LinearSegmentedColormap('WYlGn', cm_data)
WYlGn_r = LinearSegmentedColormap('WYlGn_r', reverse_cm(cm_data))

# Make a clean "from <colormap> import *"
__all__ = ['WYlGn', 'WYlGn_r']

# For use with viscm
test_cm = WYlGn

# the end
