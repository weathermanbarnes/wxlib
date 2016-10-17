#!/usr/bin/env python
# -*- encoding: utf-8

import unittest
import doctest

import dynlib.metio as m
from datetime import datetime as dt, timedelta as td


class dts2str(unittest.TestCase):
	def test_single_month_aggregation(self):
		tst = [dt(1989,1,1)+i*td(0.25) for i in range(124)]
		self.assertEqual(m.dts2str(tst, agg='six_hourly'), '198901')
		
	def test_more_than_month_aggregation(self):
		tst = [dt(1989,2,1)+i*td(0.25) for i in range(124)]
		self.assertEqual(m.dts2str(tst, agg='six_hourly'), '19890201-19890303')

	def test_short_february_aggregation(self):
		tst = [dt(1989,2,1)+i*td(0.25) for i in range(112)]
		self.assertEqual(m.dts2str(tst, agg='six_hourly'), '198902')

	def test_long_february_short_aggregation(self):
		tst = [dt(1992,2,1)+i*td(0.25) for i in range(112)]
		self.assertEqual(m.dts2str(tst, agg='six_hourly'), '19920201-19920228')

	def test_long_february_aggregation(self):
		tst = [dt(1992,2,1)+i*td(0.25) for i in range(116)]
		self.assertEqual(m.dts2str(tst, agg='six_hourly'), '199202')

	def test_december_aggregation(self):
		tst = [dt(1992,12,1)+i*td(0.25) for i in range(124)]
		self.assertEqual(m.dts2str(tst, agg='six_hourly'), '199212')

	def test_interval_mismatch(self):
		tst = [dt(1992,12,1)+i*td(0.25) for i in range(124)]
		self.assertEqual(m.dts2str(tst, agg='hourly'), '1992120100..1992123118')

	def test_uncontiguos_by_irregularity_but_right_length(self):
		tst = [dt(1992,12,1)+i*td(0.25) for i in range(124)]
		tst[17] += td(0.05)
		self.assertEqual(m.dts2str(tst, agg='six_hourly'), '199212..199212')

	def test_weekly_aggregation(self):
		tst = [dt(1992,1,1)+i*td(7) for i in range(52)]
		self.assertEqual(m.dts2str(tst, agg='weekly'), '1992')

	def test_interval_mismatch2(self):
		tst = [dt(1992,1,1)+i*td(7) for i in range(52)]
		self.assertEqual(m.dts2str(tst, agg='cal_month'), '19920101..19921223')

	def test_hourly_aggregation(self):
		tst = [dt(1986,1,1)+i*td(0,3600) for i in range(24)]
		self.assertEqual(m.dts2str(tst, agg='hourly'), '19860101')

	def test_interval_mismatch3(self):
		tst = [dt(1986,1,1)+i*td(0,3600) for i in range(24)]
		self.assertEqual(m.dts2str(tst, agg='three_hourly'), '1986010100..1986010123')

	def test_daily_two_month_aggregation(self):
		tst = [dt(1986,1,1)+i*td(1) for i in range(59)]
		self.assertEqual(m.dts2str(tst, agg='daily'), '198601-198602')

	def test_daily_multimonth_aggregation(self):
		tst = [dt(1986,1,1)+i*td(1) for i in range(181+365)]
		self.assertEqual(m.dts2str(tst, agg='daily'), '198601-198706')


# Automatically run doctests as well
def load_tests(loader, tests, ignore):
	tests.addTests(doctest.DocTestSuite(m))
	return tests



# end 
