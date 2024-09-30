#!/usr/bin/env python

''' Tendencies for meteorological variables

Currently there are routines for terms in the tendency equations of

 * Deformation and deformation angle
 * Geopotential slopes on isentropic surfaces

'''

from __future__ import absolute_import, unicode_literals

import sys

from . import dynfor
from . import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.tend, 'tend', sys.modules[__name__])



#
