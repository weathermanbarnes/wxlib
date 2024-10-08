#!/usr/bin/env python
# -*- encoding: utf-8

import os
import os.path
import subprocess
import warnings
import glob

from distutils.core import setup
from distutils.command.build_py import build_py as _build_py


# Determining dynlib version from git
version = subprocess.check_output("git describe --tags", shell=True)
version = version.strip().decode('utf8')
if "-" in version:
    idx = version.index('-')
    idx2 = version.index('-', idx+1)
    release = version[:idx]
    ncommits = version[idx+1:idx2]
    gitrev = version[idx2+1:]
    version = f'{release}.dev{ncommits}+{gitrev}'

# Are there local uncommitted changes?
changes = subprocess.check_output("git diff", shell=True)
if len(changes) > 0:
	warnings.warn('Packaging/Installing a non-committed version! Be sure you know what you are doing!')
	version += '+'

precc = 'lib/fortran/.precc'
fortran_modules = ['kind', 'config', 'consts', 'derivatives', 'detect', 'detect_lines', 'detect_rwb_contour',
		   'diag', 'ellipse', 'thermodyn', 'sphere', 'stat', 'utils']
fortran_modules = ['%s/%s.mod' % (precc, mod) for mod in fortran_modules]

# Override the build_py class to 
#  (1) make it also compile the f2py shared object 
#  (2) make python module at . the root module called dynlib
class build_py(_build_py):
	def run(self):
		subprocess.call("./compile", shell=True)
		_build_py.run(self)
		for fname in glob.glob('lib/dynfor*.so'):
			self.copy_file(fname, os.path.join(self.build_lib, 'dyn'+fname), preserve_mode=True)
		self.copy_file('lib/.dynfor_doc.pickle', os.path.join(self.build_lib, 'dynlib/.dynfor_doc.pickle'), preserve_mode=True)

		return

	def finalize_options(self):
		self.set_undefined_options('build', ('build_lib', 'build_lib'))
		_build_py.finalize_options(self)

		return

# The actual setup
setup(cmdclass={'build_py': build_py},
	name='wxlib',
	version=version,
	description='Weather systems library for gridded met data.',
	author='ARC Centre of Excellence for 21st Century Weather',
	author_email='michael.barnes@monash.edu',
	url='TBC',
	packages=['dynlib', 'dynlib.metio', 'dynlib.metio.composite', 'dynlib.plotsettings', 'dynlib.cm_dat', ],
	package_dir={'dynlib': 'lib'},
	#py_modules=['test', ],
	scripts=['bin/dynlib_init.py', ],
	data_files=[
		('lib', ['lib/libdynfor.so']),
		('include', fortran_modules),
	]
)

# the end
