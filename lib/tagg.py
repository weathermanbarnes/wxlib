#!/usr/bin/env python
# -*- encoding: utf-8

''' Time interval definitions for generating periodically averaged data

Currently the following aggregation periods exist:

 * ``all``: Temporal average everything
 * ``met_season``: Seasons after their standard meteorological definition
 * ``cal_monthly``: Months following the calendar definition
 * ``10d``: 10-day intervals starting from the given epoch date
 * ``cal_weekly``: Weeks following the calendar definition
 * ``7d`` or ``weekly``: 7-day intervals starting from the given epoch date, 
   irrespective of the day of the week
 * ``pendad``: Pentads according to their definition
 * ``5d``: 5-day intervals starting from the given epoch date
 * ``3d``: 3-day intervals starting from the given epoch date
 * ``2d``: 2-day intervals starting from the given epoch date
 * ``1d`` or ``daily``: 1-day intervals starting from the given epoch date

'''

from datetime import datetime as dt, timedelta as td
import calendar


class tagg(object):
	def __init__(self, dtstart, dtd, epoch=None, dtend=None):
		self.dtstart = dtstart
		self.dtend = dtend
		self.dtd = dtd

		if not epoch:
			self.epoch = dtstart
		else:
			self.epoch = epoch

		if dtstart == self.start(dtstart):
			self.cur = dtstart
		else: 
			self.cur = self.start_next(dtstart)
		
	def __iter__(self):
		if self.dtend:
			while self.cur < self.dtend:
				yield self.cur
				self.cur = self.start_next()
		else:
			while True:
				yield self.cur
				self.cur = self.start_next()
	
	def start(self, dti):
		raise NotImplementedError, 'To be overriden in base classes.'

	def end(self, dti):
		return self.start_next(dti) - self.dtd
	
	def start_next(self, dti=None):
		raise NotImplementedError, 'To be overriden in base classes.'




class met_season(tagg):
	def start(self, dti):
		# First month of a the season for a given month m
		fst_month = ((dti.month / 3) * 3 - 1) % 12 + 1
		# Is the beginning of the current season in a previous year?
		fst_year_offset = -1 if fst_month > dti.month else 0
		return dt(dti.year+fst_year_offset, fst_month, 1)

	def start_next(self, dti=None):
		if dti:
			dts = self.start(dti)
		else:
			dts = self.start(self.cur)
		return dt(dts.year + dts.month / 12, (dts.month + 2) % 12 + 1, 1)
	




agg = {
	'met_season': met_season,
}


# the end
