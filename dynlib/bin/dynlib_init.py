#!/usr/bin/env python
# -*- encoding: utf-8

# Step 1: Parse command line arguments to determine the name of the project directory to be created
import argparse

parser = argparse.ArgumentParser(description='Initialise a new project using dynlib')
parser.add_argument('dir', nargs=1, help='Name of the newly created project folder.')
args = parser.parse_args()
dir_name = args.dir[0]

# Step 2: Copy the default project directory from the library
import os.path
import shutil
import dynlib.dynpie

if os.path.isdir(dir_name) or os.path.islink(dir_name) or os.path.isfile(dir_name):
	raise RuntimeError, '%s already exists'

library_path = os.path.dirname(dynlib.__file__)
default_proj = os.path.join(library_path, 'default')

shutil.copytree(default_proj, dir_name)


# the end
