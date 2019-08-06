#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red': ((0.0, 1.0, 1.0), (0.33, 0.30, 0.30),  (0.867, 0.1, 0.1), (1.0, 0.5, 0.5)),
	 'green': ((0.0, 1.0, 1.0), (0.33, 0.65, 0.65),  (0.867, 0.2, 0.2), (1.0, 0.2, 0.2)),
	 'blue':  ((0.0, 1.0, 1.0), (0.33, 0.80, 0.80),  (0.867, 0.6, 0.6), (1.0, 0.8, 0.8))  }

WBuPi = LinearSegmentedColormap('WBuPi', cm_data, 256)
WBuPi_r = LinearSegmentedColormap('WBuPi_r', reverse_cm(cm_data), 256)

# Make a clean "from <colormap> import *"
__all__ = ['WBuPi', 'WBuPi_r']

# For use with viscm
test_cm = WBuPi

# the end
