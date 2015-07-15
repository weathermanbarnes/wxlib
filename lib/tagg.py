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
	''' Time interval object and iterator
	
	'''

	interval = None

	def __init__(self, dtstart, dtend=None, dtd=None, epoch=None):
		self.dtstart = dtstart
		self.dtend = dtend

		if not dtd:
			self.dtd = self.interval
		else:
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
		''' Generator function to allow iterating over the time range '''

		if self.dtend:
			while self.cur < self.dtend:
				yield self.cur
				self.cur = self.start_next()
		else:
			while True:
				yield self.cur
				self.cur = self.start_next()
	
	def start(self, dti):
		''' Start of the interval the given date is in '''

		raise NotImplementedError, 'To be overriden in base classes.'

	def end(self, dti):
		''' End date of the interval the given date is in or on '''
		
		if self.dtd <= self.interval:
			return self.start_next(dti) - self.dtd
		else:
			return self.start_next(dti)
	
	def end_after(self, dti):
		''' End date of the interval the given date is in '''
		
		if self.dtd < self.interval:
			return self.start_next(dti) - self.dtd
		else:
			return self.start_next(dti)
	
	def start_next(self, dti=None):
		''' Start of the next interval after the current or given '''

		if dti:
			dts = self.start(dti)
		else:
			dts = self.start(self.cur)

		return self._add_interval(dts)

	def start_after_or_on(self, dti):
		''' Start date of the interval starting at or after the given date '''

		dts = self.start(dti)
		if dts == dti:
			return dti
		else:
			return self._add_interval(dts)

	def _add_interval(self, dti):
		''' Add one interval to the given date '''

		return dti + self.interval




class all(tagg):
	def start(self, dti):
		__doc__ = tagg.start.__doc__

		return self.dtstart

	def _add_interval(self, dti): 
		__doc__ = tagg._add_interval.__doc__

		return dti + (self.dtend - self.dtstart)
	

class cal_year(tagg):
	interval = rtd(years=+1)

	def start(self, dti):
		__doc__ = tagg.start.__doc__

		return dt(dti.year, 1, 1)


class met_season(tagg):
	interval = rtd(months=+3)

	def start(self, dti):
		__doc__ = tagg.start.__doc__

		# First month of a the season for a given month m
		fst_month = ((dti.month / 3) * 3 - 1) % 12 + 1
		# Is the beginning of the current season in a previous year?
		fst_year_offset = -1 if fst_month > dti.month else 0
		return dt(dti.year+fst_year_offset, fst_month, 1)
	

class cal_month(tagg):
	interval = rtd(months=+1)

	def start(self, dti):
		__doc__ = tagg.start.__doc__

		return dt(dti.year, dti.month,1)


class cal_week(tagg):
	interval = rtd(weeks=+1)

	def start(self, dti):
		__doc__ = tagg.start.__doc__

		return dti + rtd(weekday=monday, hour=0, minute=0, second=0, microsecond=0)


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
		__doc__ = tagg.start.__doc__

		pentad, startleap, longpentad = self._get_pentad(dti)

		return dt(dti.year, 1, 1) + pentad*td(5) + td(startleap)

	def _add_interval(self, dti):
		__doc__ = tagg._add_interval.__doc__

		pentad, startleap, longpentad = self._get_pentad(dti)

		return dti + td(5+longpentad)


def Nday_factory(N):
	class new(tagg):
		interval = td(N)

		def start(self, dti):
			__doc__ = tagg.start.__doc__

			return self.epoch + td(((dti - self.epoch).days / N) * N)
	

	return new


def Nhour_factory(N):
	class new(tagg):
		interval = td(0, N*3600)

		def start(self, dti):
			__doc__ = tagg.start.__doc__

			return self.epoch + td(0, (int((dti - self.epoch).total_seconds()) / (N*3600)) * N*3600)
	

	return new

six_hourly = Nhour_factory(6)
daily = Nday_factory(1)
two_daily = Nday_factory(2)
three_daily = Nday_factory(3)
five_daily = Nday_factory(5)
weekly = Nday_factory(7)
ten_daily = Nday_factory(10)


__all__ = ['all', 'cal_year', 'met_season', 'cal_month', 'cal_week', 'cal_pentad', 
	'ten_daily', 'weekly', 'five_daily', 'three_daily', 'two_daily', 'daily',
	'six_hourly']

def get_by_interval(td):
	for agg in __all__:
		agg_obj = globals()[agg]
		if agg_obj.interval == td:
			return agg_obj
	
	raise ValueError, 'No aggregator found for time interval `%s`' % str(td)

__all__.append('get_by_interval')

# the end
