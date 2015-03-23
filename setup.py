#!/usr/bin/env python
# -*- encoding: utf-8

import os
import os.path
import subprocess
import warnings

from distutils.core import setup
from distutils.command.build_py import build_py as _build_py


# Determining dynlib version from git
version = subprocess.check_output("git describe --tags", shell=True)
version = version.strip()

# Are there local uncommitted changes?
changes = subprocess.check_output("git diff", shell=True)
if len(changes) > 0:
	warnings.warn('Installing a non-committed version! Be sure you know what you are doing!')
	version += '+'

# Override the build_py class to 
#  (1) make it also compile the f2py shared object 
#  (2) make python module at . the root module called dynlib
class build_py(_build_py):
	def run(self):
		subprocess.call("./compile", shell=True)
		_build_py.run(self)
		self.copy_file('dynlib.so', os.path.join(self.build_lib, 'dynlib.so'), preserve_mode=True)

		return

	def finalize_options(self):
		self.set_undefined_options('build', ('build_lib', 'build_lib'))
		self.build_lib = os.path.join(self.build_lib, 'dynlib')
		_build_py.finalize_options(self)

		return

# The actual setup
setup(cmdclass={'build_py': build_py},
	name='dynlib',
	version=version,
	description='Dynamics library for gridded met data.',
	author='Clemens Spensberger',
	author_email='csp001@uib.no',
	url='https://wikihost.uib.no/gfi/index.php/Dynlib',
	packages=['dynpie', 'tests', ],
	py_modules=['test', ],
	package_dir={'dynpie': '.'},
	package_data={'dynpie': ['default/*.py', ] },
	scripts=['scripts/dynlib_init.py', ]
)

# the end
