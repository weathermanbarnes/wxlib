#!/usr/bin/env python
# -*- encoding: utf-8

from reverse import reverse_cm
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

cm_data = {}

name = LinearSegmentedColormap(__file__, cm_data)
name_r = LinearSegmentedColormap(__file__, reverse_cm(cm_data))

# For use with viscm
test_cm = name


# the end
