#!/usr/bin/env python 
# -*- encoding: utf-8

from lib import metio as mio
import lib.context.erainterim
import lib.context.derived
import numpy as np

import unittest


class test_metsave_timeless(unittest.TestCase):
	def setUp(self):
		pass
	
	# TODO: Provide reference data to be able to split this into smaller tests
	def test_save_append_extend(self):
		s = (2,361,720)
		f, self.grid = mio.metopen('ei.ans.1979.850.T')
		f.close()
		
		# Save first
		dat = {'eqpt': np.random.random(s), 'eqpt_std': np.random.random(s) }
		try:
			mio.metsave_timeless(dat, self.grid, plev='850', name='Testme', ids=['Testme+', 'Testme-'], 
					global_atts={'cnt': 1234})
		except:
			self.fail('Failed to save initial file')
		
		# Check first
		try:
			f = mio.metopen('ei.ans.1979.Testme', mode='r+', no_static=True)
		except:
			self.fail('Failed to read initial file')
		else:
			# TODO: Add assertions to check data integrity
			f.close()
		
		# Save second (append)
		dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
		try:
			mio.metsave_timeless(dat, self.grid, plev='850', name='Testme', ids=['Testme+', 'Testme-'], 
					global_atts={'cnt': 1234})
		except:
			self.fail('Failed to save appended data to file')

		# Check second
		try:
			f = mio.metopen('ei.ans.1979.Testme', mode='r+', no_static=True)
		except:
			self.fail('Failed to read appended data')
		else:
			# TODO: Add assertions to check data integrity
			f.close()
		
		# Save second (append in new level)
		dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
		try:
			mio.metsave_timeless(dat, self.grid, plev='300', name='Testme', ids=['Testme+', 'Testme-'], 
					global_atts={'cnt': 1234})
		except:
			self.fail('Failed to save new vertial level to file')
		
		# Check second
		try:
			f = mio.metopen('ei.ans.1979.Testme', mode='r+', no_static=True)
		except:
			self.fail('Failed to read new vertical level')
		else:
			# TODO: Add assertions to check data integrity
			f.close()

		# Append should Fail: different count
		dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
		with self.assertRaises(ValueError):
			mio.metsave_timeless(dat, static, plev='500', name='Testme', ids=['Testme+', 'Testme-'], 
					global_atts={'cnt': 1235})

		# Should Fail: different shape of variables
		s_ = (1,361,720)
		dat = {'defabs': np.random.random(s_), 'defabs_min': np.random.random(s_), 'defabs_max': np.random.random(s_) }
		with self.assertRaises(ValueError):
			mio.metsave_timeless(dat, static, plev='500', name='Testme', ids=['Testme+', 'Testme-'], 
					global_atts={'cnt': 1234})

		# Should Fail: different ID names
		dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
		with self.assertRaises(ValueError):
			mio.metsave_timeless(dat, static, plev='500', name='Testme', ids=['Testme+', 'TestMe-'], 
					global_atts={'cnt': 1234})

		# Should Fail: Overwrite existing variable
		dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
		with self.assertRaises(ValueError):
			mio.metsave_timeless(dat, static, plev='300', name='Testme', ids=['Testme+', 'TestMe-'], 
					global_atts={'cnt': 1234})
		
#
