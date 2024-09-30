#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {'red':   ((0.0, 0.0, 0.0), (1.0, 0.3, 0.3)),
	   'green': ((0.0, 0.0, 0.0), (1.0, 0.6, 0.6)),
	   'blue':  ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))  }

KBu = LinearSegmentedColormap('KBu', cm_data)
KBu_r = LinearSegmentedColormap('KBu_r', reverse_cm(cm_data))

# Make a clean "from <colormap> import *"
__all__ = ['KBu', 'KBu_r']

# For use with viscm
test_cm = KBu

# the end
