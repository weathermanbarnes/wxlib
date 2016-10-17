#!/usr/bin/env python 
# -*- encoding: utf-8

import lib.metio as mio
import lib.context.erainterim
import lib.context.derived
import numpy as np


# Init
s = (2,361,720)
f, static = mio.metopen('ei.ans.1979.850.T')
f.close()

# Save first
dat = {'eqpt': np.random.random(s), 'eqpt_std': np.random.random(s) }
mio.metsave_timeless(dat, static, plev='850', name='Testme', ids=['Testme+', 'Testme-'], 
		global_atts={'cnt': 1234})

# Check first
f = mio.metopen('ei.ans.1979.Testme', mode='r+', no_static=True)
print f.dimensions.keys()
print f.variables.keys()
for group, fg in f.groups.items():
	print group, fg.variables.keys()
f.close()

# Save second (append)
dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
mio.metsave_timeless(dat, static, plev='850', name='Testme', ids=['Testme+', 'Testme-'], 
		global_atts={'cnt': 1234})

# Check second
f = mio.metopen('ei.ans.1979.Testme', mode='r+', no_static=True)
print f.dimensions.keys()
print f.variables.keys()
for group, fg in f.groups.items():
	print group, fg.variables.keys()
f.close()

# Save second (append in new level)
dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
mio.metsave_timeless(dat, static, plev='300', name='Testme', ids=['Testme+', 'Testme-'], 
		global_atts={'cnt': 1234})

# Check second
f = mio.metopen('ei.ans.1979.Testme', mode='r+', no_static=True)
print f.dimensions.keys()
print f.variables.keys()
for group, fg in f.groups.items():
	print group, fg.variables.keys()
f.close()

#
