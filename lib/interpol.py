#!/usr/bin/env python
# -*- encoding: utf-8

''' Interpolation functions '''

from __future__ import absolute_import, unicode_literals, print_function

import sys

from . import dynfor
from . import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.interpol, 'interpol', sys.modules[__name__])


#
