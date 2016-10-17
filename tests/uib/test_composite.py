#!/usr/bin/env python
# -*- encoding: utf8

from lib.composite import build, save
from lib.composite.tests import get_nh_monthly_indexes
from lib.settings import conf

import lib.context.erainterim

conf.years = range(1979,1981)

qs = [
#	('300', 'z'), 
#	('pt330', 'pv'),
	('850', 't'),
]
tests = get_nh_monthly_indexes()


stuff = build(qs, tests)
save(qs, tests, *stuff)

# the end
