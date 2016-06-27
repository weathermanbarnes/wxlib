#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import absolute_import, unicode_literals

# Make version information easily accessible
from . import dynfor

version = ''.join([x.decode('ascii') for x in dynfor.consts.version[7:]]).strip()

#
