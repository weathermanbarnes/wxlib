#!/usr/bin/env python

from __future__ import absolute_import, unicode_literals

import sys

from . import dynfor
from . import docutil

# Take over the contents of dynfor.diag_tend, dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.diag_tend, 'diag_tend', sys.modules[__name__])
docutil.takeover(dynfor.diag, 'diag', sys.modules[__name__])

#
