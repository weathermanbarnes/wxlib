#!/usr/bin/env python
# -*- encoding: utf-8

import dynlib

dynlib.diag.__doc__ = 'test'
dynlib.diag.__module__ = 'diag'
dynlib.__doc__ = 'Library of Fortran functions'

__doc__ = dynlib.__doc__
diag = dynlib.diag

# that's it
