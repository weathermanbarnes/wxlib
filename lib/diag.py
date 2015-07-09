#!/usr/bin/env python

import sys

import dynfor
import docutil

# Take over the contents of dynfor.diag_tend, dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.diag_tend, 'diag_tend', sys.modules[__name__])
docutil.takeover(dynfor.diag, 'diag', sys.modules[__name__])

#
