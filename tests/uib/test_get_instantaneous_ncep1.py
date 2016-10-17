#!/usr/bin/env python
# -*- encoding: utf8

from lib.shorthands import *
from lib.settings import conf

import lib.context.ncep1
import lib.context.plot

import unittest


class test_readplot_ncep1(unittest.TestCase):
	def test_readplot_single(self):
		try:
			dat, grid = get_instantaneous('uwnd', dt(1986,4,7,6))
		except:
			self.fail('Loading single time step NCEP data failed')
		# TODO: Add assertions to check grid
		try:
			fig.map(dat[2], grid, save=False, title=None, show=False)
		except:
			self.fail('Plotting NCEP data failed')

	def etst_read_multi(self):
		try:
			dat, grid = get_instantaneous('uwnd', (dt(1948,1,1,0), dt(1949,4,7,6)))
		except:
			self.fail('Loading single time step NCEP data failed')
		# TODO: Add assertions to check grid

#
