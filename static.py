#! /usr/bin/python
#  -*- encoding: utf-8

import numpy as np
import math

plevs  = np.array([100,200,300,400,500,550,600,650,700,750,800,850,900,950,1000])
years = range(1979,2011)

q = {'defabs': 'defabs', 'defang': 'defang',
	'u': 'u', 'v': 'v', 'Z': 'z' }
_rose = [17,]
_rose.extend(range(-18,18))
bins = { #'defabs': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22.5, 25, 27.5, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500])*1.0e-6,
	'defang': np.array(_rose)*math.pi/36.0+math.pi/72.0, }

datapath = ['/Data/gfi/work/csp001/deformation', '/scratch/reanalysis', '/media/work/reanalysis/highres', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY']
file_std   = 'ei.ans.%d.%d.%s'
file_stat  = 'ei.ans.%d.%d.%s.stat'
file_mstat = 'ei.ans.stat.%d.%s'

# Cut out    all times   all latitudes all longitudes
std_slice = (slice(None), slice(None), slice(None))

#  the end
