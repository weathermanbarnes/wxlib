#!/usr/bin/env python
# -*- encoding: utf-8

from dynlib.metio.erainterim import conf, dt, get_composite, get_static, metsave_composite
from dynlib.metio.composite import decide_by_date, decide_by_data



# Which variable(s) to composite, given by list of 2-tuples (plev, q)
qs = [('500', 'z'), ]

# Time interval to consider for the composites, here 2001-2010
timeinterval = [dt(2001,1,1,0), dt(2011,1,1,0)]


# Define the composite based on two different kinds of tests:
# 1. First define summer season
jja = decide_by_date('JJA', lambda date: date.month in [6,7,8]) 

# 2. Find time steps when Bergen 2m temperatures exceed 20 degC
warm_bgo = decide_by_data('warm_BGO', 'sfc', 't2m', lambda t2m: t2m[60,370] > 293.0)

# 3. Combine the tests
composites = [warm_bgo & jja, ]


# Actually construct the composites
dat = get_composite(qs, timeinterval, composites)

# Save the composites as netCDF files; requires injection of missing vertical level information
grid = get_static()
metsave_composite(dat, composites, grid, 'ei.ans.warm_BGO_composite')



# C'est le fin
