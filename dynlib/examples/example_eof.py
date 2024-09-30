#!/usr/bin/env python
# -*- encoding: utf8

from lib.eof import build, save, djf_mask
from lib.settings import conf

import lib.context.erainterim
import lib.context.derived

conf.years = range(1979,1989)
conf.datapath.insert(1, '/Data/gfi/spengler/csp001/jetaxis')

qs = [
#	('300', 'z'), 
#	('pt330', 'pv'),
#	('850', 't'),
	('pv2000', 'jetaxis'),
]
eofs = [
	('NH20@DJF', djf_mask, (slice(0,140), slice(None)) )
]


stuff = build(qs, eofs)
save(qs, eofs, *stuff)

# the end
