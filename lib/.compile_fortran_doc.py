#!/usr/bin/env python
# -*- encoding: utf-8

# This script ectracts documentation comments from the fortran code 
# and saves them into a pickle file. This pickle object is read in the various 
# dynlib library files in this folder (like diag.py, utils.py, etc.)
# to augment the Fortran objects on import with meaningful docstrings. 
# These docstrings are also used to create the sphinx documentation for 
# both the Fortran and python parts of dynlib.
#
# Documentation syntax: Such comments are indicated by lines beginnng with '!@'

from __future__ import unicode_literals

import sys
import docutil


fortran_doc = {}
# Expecting file names of the fortran sources as arguments
for arg in sys.argv[1:]:
	docutil.parse_fortran_file(arg, fortran_doc)

docutil.save_fortran_doc(fortran_doc)

# the end
