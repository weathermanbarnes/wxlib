#!/usr/bin/env python

import os
import sys
import pickle

import dynfor


# Take over the contents of a fortran module
pth = os.path.abspath(os.path.join(os.path.dirname(__file__), '.dynfor_doc.pickle'))
if os.path.isfile(pth):
	f = open(pth, 'r')
	fortran_doc = pickle.load(f)
	f.close()
else:
	print 'Warning: Fortran documentation not found. Make sure you have compiled dynlib.'
	fortran_doc = {}

for name, value in dynfor.ellipse.__dict__.items():
	# Backup the f2py generated docstring
	value.__f2pydoc__ = value.__doc__
	# Set the doctring from the fortran_doc extracted documentation
	value.__doc__ = fortran_doc.get('ellipse.%s' % name, None)
	# Make available as a member of this python module
	locals()[name] = value


#
