#!/usr/bin/python
# -*- encoding: utf-8

import matplotlib as mpl
mpl.use('Agg')

import use_dynlib_dev
from dynlib.shorthands import get_instantaneous, get_static, metopen, dt, plt, fig

from settings import conf, cm, proj
import composite_tests as ct

import calendar


# Should also include: EOFs
databases = {
	'testme': (dt(2011,4,7,6), dt(2011,4,25,0)),
	'testme_': dt(2012,4,7,6),
#	'testme3': (dt(2011,4,7,6), dt(2011,4,25,0), '3d'),
#	'testme3_': (dt(2012,4,7,6), '3d'),
	'testmetoo': ct.get_seasons(),
	'testmetoo_': ct.get_seasons()[0],
}

m = proj.n_hemisphere
#m = proj.N_Atlantic

plots = [
#	[('msl', 'sfc'), ],
#	[('msl_h6d', 'sfc'), ],
#	[('z', '500'), ],
#	[('z', '300'), ],
	[('t', '850'), ('z', '850'), ],
#	[('t', '500'), ('z', '500'), ('msl', 'sfc')],
#	[('pv', '300'), ],
#	[('deform', '850'), ('z', '850'), ],
#	[('deform', '500'), ('z', '500'), ],
#	[('w', 'diff_850_1000'), ],
#	[('block', 'pv2000'), ('pt', 'pv2000'), ],
#	[('rwb_a', 'pv2000'), ('pt', 'pv2000'), ],
#	[('rwb_c', 'pv2000'), ('pt', 'pv2000'), ],
#	[('rw_rfrac', '300'), ('u', '300'), ],
#	[('rw_rfrac', '500'), ('u', '500'), ],
#	[('ttr', 'sfc'), ('barbs', '850'), ],
#	[('sst', 'sfc'), ],
#	[('jetaxis', 'pv2000'), ],
#	[('jetaxis', 'pv2000'), ('pt', 'pv2000'), ],
	[('pt', 'pv2000'), ('block', 'pv2000'), ('rwb_a', 'pv2000'), ('rwb_c', 'pv2000')],
#	[('jetaxis', 'avg_900_700'), ],
]

dat = {}




#######################################################################
# Anything below should generally not be edited.

def dt2str(dt):
	return dt.strftime('%Y%m%d%H')
def dtdt2str(dta, dtz):
	if dta.hour == 0:
		if dta.day == 1:
			if dta.month == 1:
				ret = dta.strftime('%Y')
			else:
				ret = dta.strftime('%Y%m')
		else:
			ret = dta.strftime('%Y%m%d')
	else:
		ret = dta.strftime('%Y%m%d%H')
	ret += '-'

	lasthour = -conf.timestep.total_seconds()/3600 % 24
	if dta.hour == lasthour:
		if dta.day == calendar.monthrange(dta.year, dta.month)[1]:
			if dta.month == 12:
				ret += dta.strftime('%Y')
			else:
				ret += dta.strftime('%Y%m')
		else:
			ret += dta.strftime('%Y%m%d')
	else:
		ret += dta.strftime('%Y%m%d%H')

	return

def _unpack(plot):
	if len(plot) == 2:
		q, plev = plot
		plotconf = {}
	else:
		q, plev, plotconf = plot
	
	return q, plev, plotconf

def _fetch_for_composite(comp, dat, names):
	name = comp.name
	period = '%d-%d' % (conf.years[0], conf.years[-1])
	f = {}
	for plot in plots:
		for layer in plot: 
			q, plev, plotconf = _unpack(layer)
			if not (name,plev,q) in dat:
				if not plev in f:
					f[plev] = metopen(conf.file_comp % {'comp': name, 'time': period, 
							'plev': plev, 'qf': conf.qf[q]}, no_static=True)
				dat[name,plev,q] = f[plev][q+'_mean']
	
	names.append(name)

	return

static = get_static()
for database in databases:
	names = []
	# Databases is a dict
	if type(database) == str:
		conf.plotf[None,None,'name_prefix'] = database
		database = databases[database]
	
	if type(database) == dt:
		name = dt2str(database)
		for plot in plots:
			for layer in plot: 
				q, plev, plotconf = _unpack(layer)
				if (name,plev,q) not in dat:
					dat[name,plev,q], dump = get_instantaneous(q, database, plevs=plev)

		names.append(name)

	elif type(database) == ct.decider:
		_fetch_for_composite(database, dat, names)

	elif type(database) in [list, tuple]:
		if type(database[0]) == dt:
			if not type(database[-1]) == dt:
				raise NotImplementedError, '%s aggregation not available yet.' % database[-1]
			else:
				mindt = min(database)
				for plot in plots:
					for layer in plot:
						q, plev, plotconf = _unpack(layer)
						# TODO: How to check if things are already present??
						tmp, dump = get_instantaneous(q, database, plevs=plev)
						for i in range(tmp.shape[0]):
							name = dt2str(mindt + i*conf.timestep)
							dat[name,plev,q] = tmp[i,::]

		elif type(database[0]) == ct.decider:
			for comp in database:
				_fetch_for_composite(comp, dat, names)

	else:
		raise TypeError, type(database)
	
	# Loops over times, composites, eofs, etc.
	for name in names:
		# Loops over the different plots to be done for each name
		for plot in plots:
			q, plev, _plotconf = _unpack(plot[0])

			# 1. Create overlays
			overlays=[]
			if len(plot) > 1:
				for olay in plot[1:]:
					oq, oplev, oplotconf = _unpack(olay)					

					if oq == 'barbs':
						raise NotImplementedError
					elif oq == 'deform': 
						raise NotImplementedError
					elif oq == 'fronts': 
						raise NotImplementedError
					elif oq == 'jetaxis': 
						raise NotImplementedError
					else:
						overlays.append(
							fig.map_overlay_dat(dat[name,oplev,oq], static, q=oq, plev=oplev, **oplotconf)
						)

			# 2. Do the actual plot
			plotconf = conf.plotf.merge(plev, q, name=name, m=m, overlays=overlays, save='auto', show=False)
			plotconf.update(_plotconf)
			fio = fig.setup(**plotconf)
			fig.map(dat[name,plev,q], static, q=q, plev=plev, **plotconf)

			plt.close('all')

# c'est la fin!
