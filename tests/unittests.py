#!/usr/bin/env python
# -*- encoding: utf-8

import unittest
import sys
# Adding the libary path to file the modules to be tested
sys.path.append('..')

# Things to be tested
from settings import *


# Testing the settings class and instance
class settings_tester(unittest.TestCase):
	def test_object(self):
		self.assertTrue(isinstance(s, settings))
		self.assertTrue(isinstance(s.contour, settings_contour))
		
		return

	
	def test_default_settings(self):
		self.assertEqual(s.q, {'defabs': 'defabs', 'defang': 'defang', 
			'm': 'mont', 'p': 'pres', 'u': 'u', 'v': 'v', 'Z': 'z' })
		self.assertEqual(s.datapath, ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )
		self.assertEqual(s.file_std, 'ei.ans.%d.%s.%s')
		self.assertEqual(s.file_stat, 'ei.ans.%d.%s.%s.stat')
		self.assertEqual(s.file_mstat, 'ei.ans.stat.%s.%s')
		self.assertEqual(s.std_slice, slice(None))

		self.assertTrue(isinstance(s.contour, settings_contour))
		self.assertEqual(s.contour.default, {'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
			'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'ticks': None, 'ticklabels': [], 'colors': 'k', 
			'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 'levels': None, 
			'origin': None, 'extent': None, 'extend': 'neither', 'linewidths': 2.0, 'linestyles': None })
		self.assertEqual(s.contourf.v, {'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
			'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'ticks': None, 'ticklabels': [], 'colors': None, 
			'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 'levels': None, 
			'origin': None, 'extent': None, 'extend': 'neither', 'hatches': None })
		with self.assertRaises(AttributeError):
			print s.contour.deftotal

		return

	
	def test_override(self):
		s.datapath.insert(1, '/Data/gfi/scratch/csp001/deformation')
		s.datapath.insert(1, '/work/csp001/deformation')
		self.assertEqual(s.datapath, ['./', '/work/csp001/deformation', '/Data/gfi/scratch/csp001/deformation', 
			'/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )

		s.reset('datapath')
		self.assertEqual(s.datapath, ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )

		s.file_std = 'somethingelse'
		s.file_mstat = 'somethingdifferent'
		self.assertEqual(s.file_std, 'somethingelse')
		self.assertEqual(s.file_mstat, 'somethingdifferent')

		s.reset()
		self.assertEqual(s.datapath, ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )
		self.assertEqual(s.file_std, 'ei.ans.%d.%s.%s')
		self.assertEqual(s.file_stat, 'ei.ans.%d.%s.%s.stat')
		self.assertEqual(s.file_mstat, 'ei.ans.stat.%s.%s')
	
		self.assertEqual(s.contour.defabs['colors'], 'k')

		s.contour.defabs['colors'] = 'green'
		self.assertEqual(s.contour.defabs['colors'], 'green')
		self.assertEqual(s.contour.Z['colors'], 'k')

		s.contour.reset('defabs')
		self.assertEqual(s.contour.defabs['colors'], 'k')

		s.contour.Z['colors'] = 'green'
		self.assertEqual(s.contour.Z['colors'], 'green')
		s.contour.reset('Z', 'colors')
		self.assertEqual(s.contour.Z['colors'], 'k')

		return


	def test_merge(self):
		kwargs = {'colors': 'k', 'extend': 'max'}
		self.assertEqual(s.contour.merge('v', **kwargs), {'plev': None, 'lon': None, 'lat': None, 'mark': None, 
			'scale': 10, 'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'ticks': None, 'ticklabels': [],
			'colors': 'k', 'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 
			'levels': None, 'origin': None, 'extent': None, 'extend': 'max', 'linewidths': 2.0, 
			'linestyles': None })
		self.assertEqual(s.contour.v, {'plev': None, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
			'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'ticks': None, 'ticklabels': [], 'colors': 'k', 
			'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 'levels': None, 
			'origin': None, 'extent': None, 'extend': 'neither', 'linewidths': 2.0, 'linestyles': None })
		kwargs = {}
		self.assertEqual(s.contour.merge('v', **kwargs), {'plev': None, 'lon': None, 'lat': None, 'mark': None, 
			'scale': 10, 'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'ticks': None, 'ticklabels': [],
			'colors': 'k', 'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 
			'levels': None, 'origin': None, 'extent': None, 'extend': 'neither', 'linewidths': 2.0, 
			'linestyles': None })

		return
	
	
	def test_default_protection(self):
		with self.assertRaises(AttributeError):
			s.__setattr__ = object.__setattr__
		with self.assertRaises(AttributeError):
			s._settings__default = {}
		with self.assertRaises(AttributeError):
			s.contour.default = {}
		with self.assertRaises(AttributeError):
			s.contour.defabs = {}

		return



# that's it
