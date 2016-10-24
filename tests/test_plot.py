#!/usr/bin/env python
# -*- encoding: utf-8

import warnings
import unittest
import pickle
import numpy as np
import numpy.testing as nptest
from PIL import Image

from create_ref_plots import path, gfile, obfile, odfile, req_levs, ofile as reffile, conf, make_plot_list



def make_test_plot_func(name, func):
	def _tmp(self):
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			img = func(self.dat, self.grid)

		img = Image.open(img)
		imgref = Image.open(self.refplots[name])

		imgdata = np.array(img.getdata())
		imgrefdata = np.array(imgref.getdata())

		#if not np.all(imgdata == imgrefdata):
		#	img.save('test_mismatch_%s_new.png' % name)
		#	imgref.save('test_mismatch_%s_ref.png' % name)

		nptest.assert_array_equal(imgdata, imgrefdata)

		return
	
	return _tmp


class test_plot(unittest.TestCase):
	def setUp(self):
		f = open(path + '/' + gfile, 'rb')
		self.grid = pickle.load(f)
		self.grid.oro = np.zeros(self.grid.x.shape)
		self.grid.cyclic_ew = True
		f.close()
		
		self.dat = {}
		for pshort in req_levs:
			f = np.load(path + '/' + obfile % pshort)
			self.dat[pshort] = dict(f)
			f.close()
			f = np.load(path + '/' + odfile % pshort)
			self.dat[pshort].update(dict(f))
			f.close()

		f = open(path + '/' + reffile, 'rb')
		self.refplots = pickle.load(f)
		f.close()

		return


# Create test functions based on the configured diagnostics
plots = make_plot_list(save=False)
for plotname, func in plots.items():
	func_ = make_test_plot_func(plotname, func)
	setattr(test_plot, 'test_plot_'+plotname, func_)

		
#raise unittest.SkipTest('Skip this for now')
if __name__ == '__main__':
	unittest.main()

# C'est le fin
