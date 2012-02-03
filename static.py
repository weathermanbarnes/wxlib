#! /usr/bin/python
#  -*- encoding: utf-8

import numpy as np
import math

plevs  = np.array([100,200,300,400,500,550,600,650,700,750,800,850,900,950,1000])
years = range(1979,2009)

q = {'defabs': 'defabs', 'defang': 'defang',
	'u': 'u', 'v': 'v', 'Z': 'z' }
_rose = [17,]
_rose.extend(range(-18,18))
bins = {'defabs': np.array(range(50))*2.0e-5,
	'defang': np.array(_rose)*math.pi/36.0+math.pi/72.0 }

datapath = ['/Data/gfi/work/csp001/deformation', '/scratch/reanalysis', '/media/work/reanalysis/highres', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY']
file_std   = 'ei.ans.%d.%d.%s'
file_stat  = 'ei.ans.%d.%d.%s.stat'
file_mstat = 'ei.ans.stat.%d.%s'

# Cut out    all times   all latitudes all longitudes
std_slice = (slice(None), slice(None), slice(None))

#  the end
