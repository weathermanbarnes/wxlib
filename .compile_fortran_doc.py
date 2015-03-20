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


import os
import sys
import re
import pickle


# Extract the name of the subroutine/function/module from the remainder of the line
def extract_name(linefrag):
	m = re.findall('[A-Za-z0-9_]+', linefrag)
	return m[0]

# Do some plausability checks and then safe the docstring into the fortran_doc dictionary
def document(module, func, docstr):
	if not module and not func:
		raise RuntimeError, 'Got neither a module nor a function name for the docstring:\n\n%s' % docstr
	if not module:
		raise RuntimeError, 'Found a function `%s` without module' % func
	if not func:
		fortran_doc[module] = docstr
	else:
		fortran_doc['%s.%s' % (module, func)] = docstr


# Expecting file names of the fortran sources as arguments
fortran_doc = {}
for arg in sys.argv[1:]:
	if not os.path.isfile(arg):
		raise ValueError, 'File not found: `%s`' % arg

	f = open(arg)
	srclines = f.readlines()
	f.close()
	
	indoc = False
	inmodule = None
	for lnr, line in zip(range(1,len(srclines)+1), srclines):
		# Remove indent
		line = line.lstrip()

		# If line is part of a comment that is marked as part of the documentation
		if line[:2] == '!@':
			if len(line) == 2:
				_srcdoc = ''
			else:
				# Allow a little flexibility in the formatting 
				# by allowing an option whitespace character after the !@
				if line[2] in [' ', '\t']:
					_srcdoc = line[3:]
				else:
					_srcdoc = line[2:]

			if not indoc:
				indoc = True
				srcdoc = _srcdoc
			else:
				srcdoc += _srcdoc

		# If first line after a documentation comment
		elif indoc:
			known_entities = [
				'subroutine', 
				'recursive subroutine',
				'function',
				'recursive function',
			]
			documented = False
			for entity in known_entities:
				if line[:len(entity)] == entity:
					name = extract_name(line[len(entity):])
					document(inmodule, name, srcdoc)
					documented = True
			if not documented and line[:6] == 'module':
				name = extract_name(line[6:])
				document(name, None, srcdoc)
				documented = True
				inmodule = name
			if not documented:
				raise Warning, 'Found a documentation comment without associated module/function '\
						'at line %d of file `%s`, and hence will be discarded. '\
						'The docstring reads:\n\n%s' % (lnr, arg, srcdoc)
			indoc = False
			srcdoc = ''

		# Entering undocumented module 
		elif line[:6] == 'module':
			inmodule = extract_name(line[6:])


f = open('.dynlib_fortran_doc.pickle', 'w')
pickle.dump(fortran_doc, f)
f.close()

# the end
