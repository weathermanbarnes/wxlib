#!/usr/bin/env python
# -*- encoding: utf-8

from dynlib.context.erainterim import conf

from dynlib.composite import build, save
import dynlib.composite.tests as t


# Which variable(s) to composite, a given by list of tuples (plev, q)
qs = [('500', 'z'), ]

# Define the composite based on two different kinds of tests:
# 1. Get tests for the different seasons
djf, mam, jja, son = t.get_seasons()['seasons']
# 2. Define warm days in Bergen
warm_bergen = t.dat_lowerbound('warm_bergen', 't2m', 'sfc', (60, 370), 288.0)
# 3. Create composites of warm days in Bergen for spring, summer and autumn
tests = t.matrix([mam, jja, son], warm_bergen)

# Construct the composites
mean, hist, mfv, cnt, static = build(qs, tests)
# Save the composites as netCDF files.
save(qs, tests, mean, hist, mfv, cnt, static)

# C'est le fin
