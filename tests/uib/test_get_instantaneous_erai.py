#!/usr/bin/env python
# -*- encoding: utf8

from lib.shorthands import *
from lib.settings import conf

import lib.context.erainterim
import lib.context.plot

import unittest


class test_readplot_erai(unittest.TestCase):
	def test_readplot_single(self):
		try: 
			dat, grid = get_instantaneous('u', dt(1986,4,7,6), plevs=850)
		except:
			self.fail('Loading single-level/single-time step ERA-Interim data failed')
		# TODO: Add assertions to check grid
		try:
			fig.map(dat.squeeze(), grid, save=False, title=None, show=False)
		except:
			self.fail('Plotting ERA-Interim data failed')
	
	def test_read_multiple(self):
		try:
			dat, grid = get_instantaneous('u', (dt(1988,1,1,0), dt(1989,4,7,6)), plevs=['300', '500', '850'])
		except:
			self.fail('Reading multi-level/multi-time step ERA-Interim data failed')
		# TODO: Add assertions to check grid

#
