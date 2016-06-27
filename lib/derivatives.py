#!/usr/bin/env python

from __future__ import absolute_import, unicode_literals

import sys

from . import dynfor
from . import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.derivatives, 'derivatives', sys.modules[__name__])


#
