#!/usr/bin/env python
# -*- encoding: utf-8

import unittest
import pickle
import numpy as np
import numpy.testing as nptest

from create_ref_data import path, gfile, obfile, odfile, diagnostics, _prepare_args



def make_test_calc_func(pshort, resname, func, args, fourdee):
	def _tmp(self):
		dat_ = self.dat[pshort]
		# Does the diagnostic require 4d-input?
		if fourdee:
			fargs = _prepare_args(dat_, {}, self.grid, args)
			res_ = func(*fargs)
			if type(resname) == tuple:
				for resname_,res__ in zip(resname,res_):
					nptest.assert_almost_equal(res__, dat_[resname_])
			else:
				nptest.assert_almost_equal(res_, dat_[resname])
		
		# Alternatively slice by time indexes
		else:
			tlen = len(self.grid.t)
			for tidx in range(tlen):
				fargs = _prepare_args(dat_, {}, self.grid, args, slc=(tidx, slice(None), slice(None), slice(None)))
				res_ = func(*fargs)
				if type(resname) == tuple:
					for resname_,res__ in zip(resname,res_):
						nptest.assert_almost_equal(res__, dat_[resname_][tidx,:,:,:])
				else:
					nptest.assert_almost_equal(res_, dat_[resname][tidx,:,:,:])

		return
	
	return _tmp


class test_func(unittest.TestCase):
	def setUp(self):
		f = open(path + '/' + gfile, 'rb')
		self.grid = pickle.load(f)
		f.close()
		
		self.dat = {}
		for pshort in diagnostics:
			f = np.load(path + '/' + obfile % pshort)
			self.dat[pshort] = dict(f)
			f.close()
			f = np.load(path + '/' + odfile % pshort)
			self.dat[pshort].update(dict(f))
			f.close()

		return


# Create test functions based on the configured diagnostics
for pshort, diags in diagnostics.items():
	for resname, func, args, fourdee in diags:
		func_ = make_test_calc_func(pshort, resname, func, args, fourdee)

		if type(resname) == tuple:
			resname = resname[0]
		setattr(test_func, 'test_calc_'+resname, func_)

		
#raise unittest.SkipTest('Skip this for now')


# C'est le fin
