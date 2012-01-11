#! /usr/bin/python
#  -*- encoding: utf-8

plevs  = [100,200,300,400,500,550,600,650,700,750,800,850,900,950,1000]
years = range(1979,1989)

q = {'defabs': 'defabs', 'defang': 'defang',
	'u': 'u', 'v': 'v', 'Z': 'z' }

datapath = ['/scratch/reanalysis', '/media/work/reanalysis/highres', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY']
file_std   = 'ei.ans.%d.%d.%s'
file_stat  = 'ei.ans.%d.%d.%s.stat'
file_mstat = 'ei.ans.stat.%d.%s'

# Cut out    all times    lat > 15°N       all longitudes
std_slice = (slice(None), slice(None,151), slice(None))

#  the end
