#!/usr/bin/env python
# -*- encoding: utf-8

from dynlib.dynpie.settings import conf as c
import dynlib.dynpie.proj as proj
import dynlib.dynpie.cm as cm


# Adapt the default settings to your needs here, then use 
# >>> from settings import *
# in your personal scripts. 
#
# Examples:
# (1) Add a folder to the datapath which is scanned in metopen. 
#     The first argument is the location within the path, the second place in the example below.
# >>> c.datapath.insert(1, '/work/csp001/deformation')
# (2) Change the default colormap used in contourf-plots for the quantity 'Z'
# >>> c.contourf.Z['cmap'] = plt.cm.RdBu_r
#
# Note that you can always access the default settings by
# >>> from dynlib import default_conf
# >>> default_conf
# >>> print default_conf.datapath
# >>> print default_conf.contourf.Z

# that's it
