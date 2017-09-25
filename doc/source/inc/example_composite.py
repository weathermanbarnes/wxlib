#!/usr/bin/env python
# -*- encoding: utf-8

from dynlib.context.erainterim import conf

from dynlib.composite import build, save
import dynlib.composite.tests as t


# Which variable(s) to composite, given by list of 2-tuples (plev, q)
qs = [('500', 'z'), ]

# Define the composite based on two different kinds of tests:
# 1. First, get tests for the different seasons
djf, mam, jja, son = t.get_seasons()['seasons']
# 2. Then tests based on the monthly values common climate indices 
#    for the northern hemisphere, e.g. NAO and PNA
tests = t.get_nh_monthly_indexes()
# 3. Combine the tests to make composites for all combinations between
#    NAO+/-, PNA+/- for both summer and winter
tests = t.matrix([djf, jja,], tests)

# Actually construct the composites
mean, hist, mfv, cnt, static = build(qs, tests)
# Save the composites as netCDF files.
save(qs, tests, mean, hist, mfv, cnt, static)

# C'est le fin
