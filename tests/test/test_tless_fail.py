#!/usr/bin/env python 
# -*- encoding: utf-8

import lib.metio as mio
import lib.context.erainterim
import lib.context.derived
import numpy as np
from copy import copy

# Init
s = (2,361,720)
f, static = mio.metopen('ei.ans.1979.850.T')
f.close()


if False:
	# Should Fail: different count
	dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
	mio.metsave_timeless(dat, static, plev='500', name='Testme', ids=['Testme+', 'Testme-'], 
			global_atts={'cnt': 1235})

if False:
	# Should Fail: different shape of variables
	s_ = (1,361,720)
	dat = {'defabs': np.random.random(s_), 'defabs_min': np.random.random(s_), 'defabs_max': np.random.random(s_) }
	mio.metsave_timeless(dat, static, plev='500', name='Testme', ids=['Testme+', 'Testme-'], 
			global_atts={'cnt': 1234})

if False:
	# Should Fail: different ID names
	dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
	mio.metsave_timeless(dat, static, plev='500', name='Testme', ids=['Testme+', 'TestMe-'], 
			global_atts={'cnt': 1234})

if False:
	# Should Fail: Overwrite existing variable
	dat = {'defabs': np.random.random(s), 'defabs_min': np.random.random(s), 'defabs_max': np.random.random(s) }
	mio.metsave_timeless(dat, static, plev='300', name='Testme', ids=['Testme+', 'TestMe-'], 
			global_atts={'cnt': 1234})

#
