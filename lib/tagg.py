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
from dateutil.relativedelta import relativedelta as rtd, MO as monday
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
		if dti:
			dts = self.start(dti)
		else:
			dts = self.start(self.cur)

		return self._add_interval(dts)

	def _add_interval(self, dti):
		raise NotImplementedError, 'To be overriden in base classes.'




class all(tagg):
	def start(self, dti):
		return self.dtstart

	def _add_interval(self, dti): 
		return dti + (self.dtend - self.dtstart)
	

class met_season(tagg):
	def start(self, dti):
		# First month of a the season for a given month m
		fst_month = ((dti.month / 3) * 3 - 1) % 12 + 1
		# Is the beginning of the current season in a previous year?
		fst_year_offset = -1 if fst_month > dti.month else 0
		return dt(dti.year+fst_year_offset, fst_month, 1)

	def _add_interval(self, dti): 
		return dti + rtd(months=+3)
	

class cal_month(tagg):
	def start(self, dti):
		return dt(dti.year, dti.month,1)

	def _add_interval(self, dti):
		return dti + rtd(months=+1)


class cal_week(tagg):
	def start(self, dti):
		return dti + rtd(weekday=monday, hour=0, minute=0, second=0, microsecond=0)

	def _add_interval(self, dti):
		return dti + rtd(weeks=+1)


class cal_pentad(tagg):
	def _get_pentad(self, dti):
		''' Calculate pentad number (start counting with 0), 
		if pentad start date is after Feb 29 on a leap year,
		and if the the pentad includes the leap day
		'''

		startleap = 0
		longpentad = 0
		if calendar.isleap(dti.year) and dti.timetuple().tm_yday >= 60:
			startleap = 1

		pentad = ((dti - dt(dti.year, 1, 1)).days - startleap)/ 5 
		if calendar.isleap(dti.year) and pentad == 11:
			longpentad = 1

		return pentad, startleap, longpentad

	def start(self, dti):
		pentad, startleap, longpentad = self._get_pentad(dti)

		return dt(dti.year, 1, 1) + pentad*td(5) + td(startleap)

	def _add_interval(self, dti):
		pentad, startleap, longpentad = self._get_pentad(dti)

		return dti + td(5+longpentad)


def Nday_factory(N):
	class new(tagg):
		def start(self, dti):
			return self.epoch + td(((dti - self.epoch).days / N) * N)

		def _add_interval(self, dti):
			return dti + td(N)
	
	return new

daily = Nday_factory(1)
two_daily = Nday_factory(2)
three_daily = Nday_factory(3)
five_daily = Nday_factory(5)
weekly = Nday_factory(7)
ten_daily = Nday_factory(10)


__all__ = ['met_season', 'cal_month', 'cal_week', 'cal_pentad', 
	'ten_daily', 'weekly', 'five_daily', 'three_daily', 'two_daily', 'daily']


# the end
