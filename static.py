#! /usr/bin/python
#  -*- encoding: utf-8

plevs  = ['pt350'] #['pt350','pt330','pt315','pt300']
years = [1981]
#years = range(1979,2011) #nb end with final year+1
#plevs=['pt350']
#years=[1979]

q = {'defabs': 'defabs', 'defang': 'defang',
	'u': 'u', 'v': 'v', 'Z': 'z' }

datapath = ['/Data/gfi/users/tsp065/students/amu006/clim', '/Data/gfi/users/tsp065/students/amu006/rwbdata', '/Data/gfi/work/csp001/deformation', '/scratch/reanalysis', '/media/work/reanalysis/highres', '/Data/gfi/share/Reanalyses/ERA_INTERIM/6HOURLY']
file_std   = 'ei.ans.%d.%d.%s'
file_stat  = 'ei.ans.%d.%d.%s.stat'
file_mstat = 'ei.ans.stat.%d.%s'

file_clim = 'ei.ans.%s.%s_clim' # eg (plev, 'strain')

# Cut out     times        lats         longitudes
std_slice = (slice(None), slice(None), slice(None))

#  the end

# Useful Slices:
# start 29 jan 0600  slice(114,313)
# jan and feb slice(None,232)
