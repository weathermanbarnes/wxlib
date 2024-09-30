#!/usr/bin/env python
# -*- encoding: utf-8

from dynlib.metio.erainterim import conf, dt, get_composite, get_static, metopen, metsave_composite
from dynlib.metio.composite import decide_by_date, decide_by_timeseries, matrix



# Which variable(s) to composite, given by list of 2-tuples (plev, q)
qs = [('500', 'z'), ]

# Time interval to consider for the composites, here 2001-2010
timeinterval = [dt(2001,1,1,0), dt(2011,1,1,0)]


# Define the composite based on two different kinds of tests:
# 1. First, get tests for the different seasons
djf = decide_by_date('DJF', lambda date: date.month in [12,1,2])
jja = decide_by_date('JJA', lambda date: date.month in [6,7,8]) 

# 2. Then tests based on the monthly values of the NAO
ts = metopen('indexes/ts_NAO', no_static=True)
naop = decide_by_timeseries('NAO+', ts, lambda idx: idx >= 1.0)
naom = decide_by_timeseries('NAO-', ts, lambda idx: idx <= -1.0)

# 3. Combine the tests to make composites for all combinations between
#    NAO+/- and summer/winter
composites = matrix([djf, jja,], [naom, naop])


# Actually construct the composites
dat = get_composite(qs, timeinterval, composites)

# Save the composites as netCDF files; requires injection of missing vertical level information
grid = get_static()
metsave_composite(dat, composites, grid, 'ei.ans.NAO_composites')



# C'est le fin
