#!/usr/bin/env python

import sys

import dynfor
import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.derivatives, 'derivatives', sys.modules[__name__])


#
