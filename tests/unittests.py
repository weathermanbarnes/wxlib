#!/usr/bin/env python
# -*- encoding: utf-8

import unittest
import sys

# Things to be tested
from dynpie.settings import *


# Testing the settings class and instance
class settings_tester(unittest.TestCase):
	def test_object(self):
		self.assertTrue(isinstance(conf, settings))
		self.assertTrue(isinstance(conf.contour, settings_contour))
		
		return

	
	def test_default_settings(self):
		self.maxDiff = None
		self.assertEqual(conf.q, {'defabs': 'defabs', 'defang': 'defang', 
			'm': 'mont', 'p': 'pres', 'u': 'u', 'v': 'v', 'Z': 'z', 'oro': 'oro' })
		self.assertEqual(conf.datapath, ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )
		self.assertEqual(conf.file_std, 'ei.ans.%d.%s.%s')
		self.assertEqual(conf.file_stat, 'ei.ans.%d.%s.%s.stat')
		self.assertEqual(conf.file_mstat, 'ei.ans.stat.%s.%s')
		self.assertEqual(conf.std_slice, (slice(None), slice(None), slice(None)))

		self.assertTrue(isinstance(conf.contour, settings_contour))
		self.assertEqual(conf.contour.default, {'m': wmap, 'plev': 800, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
			'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'oroscale': scale_oro, 'ticks': None, 'ticklabels': [], 'colors': 'k', 
			'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 'levels': None, 
			'origin': None, 'extent': None, 'extend': 'neither', 'linewidths': 2.0, 'linestyles': None })
		self.assertEqual(conf.contourf.v, {'m': wmap, 'plev': 800, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
			'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'oroscale': scale_oro, 'ticks': None, 'ticklabels': [], 'colors': None, 
			'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 'levels': None, 
			'origin': None, 'extent': None, 'extend': 'neither', 'hatches': None })
		with self.assertRaises(AttributeError):
			print conf.contour.deftotal

		return

	
	def test_override(self):
		conf.datapath.insert(1, '/Data/gfi/scratch/csp001/deformation')
		conf.datapath.insert(1, '/work/csp001/deformation')
		self.assertEqual(conf.datapath, ['./', '/work/csp001/deformation', '/Data/gfi/scratch/csp001/deformation', 
			'/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )

		conf.reset('datapath')
		self.assertEqual(conf.datapath, ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )

		conf.file_std = 'somethingelse'
		conf.file_mstat = 'somethingdifferent'
		self.assertEqual(conf.file_std, 'somethingelse')
		self.assertEqual(conf.file_mstat, 'somethingdifferent')

		conf.reset()
		self.assertEqual(conf.datapath, ['./', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY'] )
		self.assertEqual(conf.file_std, 'ei.ans.%d.%s.%s')
		self.assertEqual(conf.file_stat, 'ei.ans.%d.%s.%s.stat')
		self.assertEqual(conf.file_mstat, 'ei.ans.stat.%s.%s')
	
		self.assertEqual(conf.contour.defabs['colors'], 'k')

		conf.contour.defabs['colors'] = 'green'
		self.assertEqual(conf.contour.defabs['colors'], 'green')
		self.assertEqual(conf.contour.Z['colors'], 'k')

		conf.contour.reset('defabs')
		self.assertEqual(conf.contour.defabs['colors'], 'k')

		conf.contour.Z['colors'] = 'green'
		self.assertEqual(conf.contour.Z['colors'], 'green')
		conf.contour.reset('Z', 'colors')
		self.assertEqual(conf.contour.Z['colors'], 'k')

		return


	def test_merge(self):
		self.maxDiff = None
		kwargs = {'colors': 'k', 'extend': 'max'}
		self.assertEqual(conf.contour.merge('v', **kwargs), {'m': wmap, 'plev': 800, 'lon': None, 'lat': None, 'mark': None, 
			'scale': 10, 'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'oroscale': scale_oro, 'ticks': None, 'ticklabels': [],
			'colors': 'k', 'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 
			'levels': None, 'origin': None, 'extent': None, 'extend': 'max', 'linewidths': 2.0, 
			'linestyles': None })
		self.assertEqual(conf.contour.v, {'m': wmap, 'plev': 800, 'lon': None, 'lat': None, 'mark': None, 'scale': 10, 
			'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'oroscale': scale_oro, 'ticks': None, 'ticklabels': [], 'colors': 'k', 
			'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 'levels': None, 
			'origin': None, 'extent': None, 'extend': 'neither', 'linewidths': 2.0, 'linestyles': None })
		kwargs = {}
		self.assertEqual(conf.contour.merge('v', **kwargs), {'m': wmap, 'plev': 800, 'lon': None, 'lat': None, 'mark': None, 
			'scale': 10, 'overlays': [], 'disable_cb': False, 'show': True, 'save': '', 'title': '', 
			'coastcolor': 'k', 'gridcolor': 'k', 'maskcolor': '0.25', 'orocolor': 'k', 'oroalpha': 0.4,
			'oroscale': scale_oro, 'ticks': None, 'ticklabels': [],
			'colors': 'k', 'alpha': 1.0, 'cmap': None, 'norm': None, 'vmin': None, 'vmax': None, 
			'levels': None, 'origin': None, 'extent': None, 'extend': 'neither', 'linewidths': 2.0, 
			'linestyles': None })

		return
	
	
	def test_default_protection(self):
		with self.assertRaises(AttributeError):
			conf.__setattr__ = object.__setattr__
		with self.assertRaises(AttributeError):
			conf._settings__default = {}
		with self.assertRaises(AttributeError):
			conf.contour.default = {}
		with self.assertRaises(AttributeError):
			conf.contour.defabs = {}

		return



# that's it
