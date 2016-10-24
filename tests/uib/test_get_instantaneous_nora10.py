#!/usr/bin/env python
# -*- encoding: utf8

from lib.shorthands import *
from lib.settings import conf
from lib import proj

import lib.context.nora10
import lib.context.plot

import unittest


class test_readplot_nora10(unittest.TestCase):
	def test_readplot_single(self):
		try: 
			dat, grid = get_instantaneous('u', dt(1986,4,7,6), plevs=850)
		except:
			self.fail('Loading single-level/single-time step NORA10 data failed')
		# TODO: Add assertions to check grid
		try:
			fig.map(dat.squeeze(), grid, save=False, show=False, title=None, m=proj.N_Atlantic)
		except:
			self.fail('Plotting NORA10 data failed')
	
	def test_read_multi(self):
		try:
			dat, grid = get_instantaneous('u', (dt(1988,1,1,0), dt(1989,4,7,6)), plevs=['300', '500', '850'])
		except:
			self.fail('Loading multi-level/multi-time step NORA10 data failed')
		# TODO: Add assertions to check grid

#
