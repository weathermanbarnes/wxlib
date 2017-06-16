#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import absolute_import, unicode_literals

# Make version information easily accessible
from . import dynfor

# Different version of numpy/f2py behave differently: 
# Version string is either array of character or scalar string
if len(dynfor.consts.version.shape) == 1:
	version = ''.join([x.decode('ascii') for x in dynfor.consts.version[7:]]).strip()
else:
	version = str(dynfor.consts.version)[9:-1].strip()


#
