#! /usr/bin/python
#  -*- encoding: utf-8

plevs  = ['pt315']
years = [2009]

q = {'defabs': 'defabs', 'defang': 'defang',
	'u': 'u', 'v': 'v', 'Z': 'z' }

datapath = ['/scratch/reanalysis', '/media/work/reanalysis/highres', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY']
file_std   = 'ei.ans.%d.%d.%s'
file_stat  = 'ei.ans.%d.%d.%s.stat'
file_mstat = 'ei.ans.stat.%d.%s'

# Cut out    start feb 09   all lats       all longitudes
std_slice = (slice(125,145), slice(None), slice(None))

#  the end
