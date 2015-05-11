#!/usr/bin/env python

''' Statistics functions

Collection of statistical utilities meant to be applied to a
(potentially larger) gridded data set.
'''

import sys

import dynfor
import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.stat, 'stat', sys.modules[__name__])





#
