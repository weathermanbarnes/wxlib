#! /usr/bin/python

plevs  = [100,200,300,400,500,550,600,650,700,750,800,850,900,950,1000]
years = range(1979,1989)

q = {'defabs': 'da', 'defang': 'da', 'deform': 'd', 'defnor': 'dn',
	'u': 'u', 'v': 'v', 'Z': 'Z' }

datapath = ['/scratch/reanalysis', '/media/work/reanalysis/highres']
file_std   = 'ei.ans.%d.%d.%s.npz'
file_stat  = 'ei.ans.%d.%d.%s.stat.npz'
file_mstat = 'ei.ans.stat.%d.%s.npz'

#  the end
