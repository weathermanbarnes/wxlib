#!/usr/bin/env python
# -*- encoding: utf8


from dynlib.settings import conf
from dynlib.metio import metopen

__metafun__ = {}

class dynfun(object):
	def __init__(self, fun, dependency_tables, provides_tables, uses_dx=False, uses_dy=False, uses_dz=False):
		self.deps = dependency_tables
		self.provides = provides_tables
		self.fun = fun
		
		# Register what this function can do, to be able to recursively resolve dependencies
		for table in provides_tables:
			for plev in table[0]:
				for q in table[1]: 
					if (plev, q) in __metafun__:
						__metafun__[plev,q].append(self)
					else:
						__metafun__[plev,q] = [self, ]

		self.uses_dx = uses_dx
		self.uses_dy = uses_dy
		self.uses_dz = uses_dz
		
		return


	def __call__(self, plev_in=None):
		# Dependencies might be specified using plev=None as a wild card if applicable on all 
		# vertical levels. In this case, plev_in must be specified such that the function knows
		# what to calculate.
		specific_deps = []
		for table in self.deps:
			specific_plevs = []
			for plev in table[0]:
				if plev == None:
					if plev_in == None:
						raise ValueError, 'This dynlib function needs to know which vertical level to operate on, but None was given.'
					if not plev_in in specific_plevs:
						specific_plevs.append(plev_in)

			specific_deps.append((specific_plevs, table[1]))
		
		for year in conf.years:
			# Fetch data
			args = []
			for table in specific_deps:
				for plev in table[0]:
					for q in table[1]:
						f, dat, static = metopen(conf.file_std % {
							'time': year, 'plev': plev, 'qf': conf.qf[q]}, q)
						if dat in LINES:
							lq = LINES[q]
							f, datoff = metopen(conf.file_std % {
								'time': year, 'plev': plev, 'qf': conf.qf[lq]}, lq)
							dat = dynlib.utils.mask_lines(dat, datoff)

						args.append(dat)

			if self.uses_dx:
				args.append(static.dx)
			if self.uses_dy:
				args.append(static.dy)
			if self.uses_dz:
				args.append(static.dz)
			
			# Calculate and yield result
			yield self.fun(*args)
		


# TODO: Status quo: 
#
# Need to build a database of what each function uses and provides, in order to
# (a) Define what commands exist using dyncal
# (b) Know which dynlib function to call which what arguments
# (c) Potentially build dependency trees for chains of functions to call for a specific task
#
# Furthermore, a few additional commands must be available: e.g. "plot". These additional commands
# cannot be directly mapped to a dynlib function (or maybe they should be?), but to supplementary 
# scripts like dyncal_plot.py

# the end
