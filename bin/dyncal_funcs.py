#!/usr/bin/env python
# -*- encoding: utf8


from dynlib.settings import conf
from dynlib.metio import metopen

import dynlib.diag


# TODO: LINES should be centralised somewhere in the variable definitions!
LINES = {'fronts': 'froff', 'convl': 'cloff', 'defl': 'dloff', 'vorl': 'vloff', 'jetaxis': 'jaoff'}


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


	def __call__(self, plev_in=None, q_in=None):
		# Dependencies might be specified using plev=None as a wild card if applicable on all 
		# vertical levels. In this case, plev_in must be specified such that the function knows
		# what to calculate.

		q_fill = ()

		specific_deps = []
		for table in self.deps:
			specific_plevs = []
			for plev in table[0]:
				if plev == None:
					if plev_in == None:
						raise ValueError, 'This dynlib function needs to know which vertical level to operate on, but None was given.'
					if not plev_in in specific_plevs:
						specific_plevs.append(plev_in)
				else: 
					specific_plevs.append(plev)

			specific_q = []
			for q in table[1]:
				if q == None: 
					if q_in == None:
						raise ValueError, 'This dynlib function needs to know with data to operate on, but None was given.'
					q_fill = (q_in, )
					if not q_in in specific_q:
						specific_q.append(q_in)
				else:
					specific_q.append(q)

			specific_deps.append((specific_plevs, specific_q))
		
		# Fill in blanks in the provides section for functions with variable input data
		if q_fill:
			self.specific_provides = []
			for table in self.provides:
				specific_q = []
				for q in table[1]:
					if '%s' in q:
						specific_q.append(q % q_fill)
					else:
						specific_q.append(q)
				self.specific_provides.append((table[0], specific_q))
		else:
			self.specific_provides = self.provides
		
		for year in conf.years:
			# Fetch data
			args = []
			for table in specific_deps:
				for plev in table[0]:
					for q in table[1]:
						f, dat, static = metopen(conf.file_std % {
							'time': year, 'plev': plev, 'qf': conf.qf[q]}, q)
						if q in LINES:
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
			yield self.fun(*args), static


	
# Vorticity and divergence
cal_vo = dynfun(dynlib.diag.vor, [([None,], ['u', 'v',]), ], [([None,], ['vo',]), ], uses_dx=True, uses_dy=True)
cal_div = dynfun(dynlib.diag.div, [([None,], ['u', 'v',]), ], [([None,], ['div',]), ], uses_dx=True, uses_dy=True)

# Deformation
cal_def_st = dynfun(dynlib.diag.def_stretch, [([None,], ['u', 'v',]), ], [([None,], ['def_st',]), ], uses_dx=True, uses_dy=True)
cal_def_sh = dynfun(dynlib.diag.def_shear, [([None,], ['u', 'v',]), ], [([None,], ['def_sh',]), ], uses_dx=True, uses_dy=True)
cal_defabs = dynfun(dynlib.diag.def_total, [([None,], ['u', 'v',]), ], [([None,], ['defabs',]), ], uses_dx=True, uses_dy=True)
cal_defang = dynfun(dynlib.diag.def_angle, [([None,], ['u', 'v',]), ], [([None,], ['defang',]), ], uses_dx=True, uses_dy=True)
cal_defanr = dynfun(dynlib.diag.def_angle_nat, [([None,], ['u', 'v',]), ], [([None,], ['defanr',]), ], uses_dx=True, uses_dy=True)
cal_def_natshl = dynfun(dynlib.diag.def_nat_shearless, [([None,], ['u', 'v',]), ], [([None,], ['defabs_shl', 'defanr_shl',]), ], uses_dx=True, uses_dy=True)
cal_ow = dynfun(dynlib.diag.okuboweiss, [([None,], ['u', 'v',]), ], [([None,], ['ow',]), ], uses_dx=True, uses_dy=True)

# Frontogenesis
cal_frontogenesis = dynfun(dynlib.diag.frontogenesis, [([None,], ['u', 'v', None]), ], [([None,], ['%s_stretch', '%s_stir',]), ], uses_dx=True, uses_dy=True)


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
