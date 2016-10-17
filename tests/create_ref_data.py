#!/usr/bin/env python3

import numpy as np
import pygrib
import pickle
import os.path

from dynlib import dynfor
import dynlib.derivatives
import dynlib.detect
import dynlib.diag
import dynlib.humidity
import dynlib.tend


Rl = dynfor.consts.rl


''' Read CFSR grib files and extract relevant base data for dynlib testing '''


def calc(dat, res, grid, diagnostics, slicebylev=True):
	''' Generic diagnostics calculator 
	
	This function prepares the arguments to be given to each calculation routine
	and takes care of the dimensionality (3D vs 4D) of the expected arguments.
	'''

	# Prepare arguments
	def _prepare_args(args, slc=slice(None)):
		fargs = []
		for var in args:
			if var == 'dx':
				fargs.append(grid.dx)
			elif var == 'dy':
				fargs.append(grid.dy)
			elif var == 'lvl':
				fargs.append(dat[var][slc[1]])
			elif var == 'fcor':
				if grid.gridtype == 'latlon':
					fargs.append(4*np.pi*np.cos(grid.y*np.pi/180.0)/86164)
				else:
					fargs.append(dat[var][slc[2:]])
			elif var == 'lat':
				if grid.gridtype == 'latlon':
					fargs.append(grid.y[:,0])
				else:
					fargs.append(dat[var][slc[2:3],0])
			elif type(var) not in [str, ]:
				fargs.append(var)
			elif var in dat:
				fargs.append(dat[var][slc])
			elif var in res:
				fargs.append(res[var][slc])
			else:
				raise ValueError('Unknown variable %s' % var)

		return fargs
	
	# Do the actual calculations
	for resname, func, args, fourdee in diagnostics:
		print('Calculating %s' % str(resname))

		# Does the diagnostic require 4d-input?
		if fourdee:
			fargs = _prepare_args(args)
			res_ = func(*fargs)
			if type(resname) == tuple:
				for resname_,res__ in zip(resname,res_):
					res[resname_] = res__
			else:
				res[resname] = res_
		
		# if not, slice the input into vertical slices
		elif slicebylev:
			if type(resname) == tuple:
				for resname_ in resname:
					res[resname_] = None
			else:
				res[resname] = None
			
			plen = len(grid.z)
			for pidx in range(plen):
				fargs = _prepare_args(args, slc=(slice(None), pidx, slice(None), slice(None)))
				res_ = func(*fargs)
				if type(resname) == tuple:
					for resname_,res__ in zip(resname,res_):
						if type(res[resname_]) == type(None):
							s_ = res_.shape
							res[resname_] = np.empty(s_[0:1] + (plen,) + s_[1:])
						res[resname_][:,pidx,:,:] = res__
				else:
					if type(res[resname]) == type(None):
						s_ = res_.shape
						res[resname] = np.empty(s_[0:1] + (plen,) + s_[1:])
					res[resname][:,pidx,:,:] = res_

		# Alternatively slice by time indexes
		else:
			if type(resname) == tuple:
				for resname_ in resname:
					res[resname_] = None
			else:
				res[resname] = None
			
			tlen = len(grid.t)
			for tidx in range(tlen):
				fargs = _prepare_args(args, slc=(tidx, slice(None), slice(None), slice(None)))
				res_ = func(*fargs)
				if type(resname) == tuple:
					for resname_,res__ in zip(resname,res_):
						if type(res[resname_]) == type(None):
							s_ = res__.shape
							res[resname_] = np.empty((tlen,) + s_)
						res[resname_][tidx,:,:,:] = res__
				else:
					if type(res[resname]) == type(None):
						s_ = res_.shape
						res[resname] = np.empty((tlen,) + s_)
					res[resname][tidx,:,:,:] = res_

	return res


def make_pressure(plevs, anydat):
	''' Create a 4d pressure field from pressure level information '''
	return np.ones(anydat.shape) * np.array(plevs)[:,np.newaxis,np.newaxis]

def density(t, p):
	''' Calculate density from ideal gas law '''
	return Rl*t/p

def dz(dat):
	datz = np.empty(dat.shape)*np.nan
	datz[1:-1,:,:] = dat[2:,:,:] - dat[:-2,:,:]
	return datz

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ifile = 'cdas1.t00z.pgrbh00.grib2'
gfile = 'cfsr_grid.pickle'
obfile = 'ref_data_base_%s.npz'
odfile = 'ref_data_derived_%s.npz'


# Keys and values in this dict correspond to the grib attributes
# - nameOfFirstFixedSurface (vertical level name)
# - scaledValueOfFirstFixedSurface (vertical level value)
# - shortNameECMF (sic! short variable name follwing ECMWF convention)
base_data = {
	'Isobaric surface': ('p', [82500, 85000, 87500, ], 
		['u', 'v', 'w', 'gh', 't', 'q',	], 
	),
	'Potential vorticity surface': ('pv', [2000, ], 
		['u', 'v', 't', 'pres', ], 
	),
}

# Currently untested: 
# - dynlib.diag.geop_from_montgp (missing input data)
diagnostics = { 'pv': [
		('defang', dynlib.diag.def_angle, ['u', 'v', 'dx', 'dy'], False),
		('defang_nat', dynlib.diag.def_angle_nat, ['u', 'v', 'dx', 'dy'], False),
		(('defabs_nat_shearless', 'defang_nat_shearless'), 
			dynlib.diag.def_nat_shearless, ['u', 'v', 'dx', 'dy'], False),
		('def_shear', dynlib.diag.def_shear, ['u', 'v', 'dx', 'dy'], False),
		('def_stretch', dynlib.diag.def_stretch, ['u', 'v', 'dx', 'dy'], False),
		('defabs', dynlib.diag.def_total, ['u', 'v', 'dx', 'dy'], False),
		('div', dynlib.diag.div, ['u', 'v', 'dx', 'dy'], False),
		('pt', dynlib.humidity.theta, ['t', 'pres'], False),
		(('fge', 'fstir'), dynlib.diag.frontogenesis, ['pt', 'u', 'v', 'dx', 'dy'], False),
		('isoline_angle', dynlib.diag.isoline_angle, ['pt', 'dx', 'dy'], False),
		('isoline_to_deformation_angle', dynlib.diag.isoline_to_deformation_angle, ['u', 'v', 'pt', 'dx', 'dy'], False),
		('okuboweiss', dynlib.diag.okuboweiss, ['u', 'v', 'dx', 'dy'], False),
		('shear_nat', dynlib.diag.shear_nat, ['u', 'v', 'dx', 'dy'], False),
		('vor', dynlib.diag.vor, ['u', 'v', 'dx', 'dy'], False),
	], 
	'p': [
		('pres', make_pressure, ['lvl', 't'], False),
		('dp', dz, ['pres'], False),
		('rho', density, ['t', 'pres'], False),
		(('accgrad_eigpr', 'accgrad_eigpi', 'accgrad_eigmr', 'accgrad_eigmi'), 
			dynlib.diag.accgrad_eigs, ['u', 'v', 'gh', 'lat', 'dx', 'dy'], False),
		('rsr', dynlib.diag.rotation_strain_ratio, ['u', 'v', 'gh', 'fcor', 'dx', 'dy'], False),
		(('ug', 'vg'), dynlib.diag.uv_geo_from_pot, ['gh', 'lat', 'dx', 'dy'], False),
		(('defabs3d', 'defaxes3d'), dynlib.diag.def_3d, ['u', 'v', 'w', 'rho', 'dx', 'dy', 'dp', False, False], True),
	],
}



if __name__ == '__main__':
	bdat = {}	# Base data (read from grib file)
	ddat = {}	# Derived data (reference calculations by dynlib)

	f = open(gfile, 'rb')
	grid = pickle.load(f)
	f.close()

	# Read file: Either cached base data, or from grib file
	f = pygrib.open(ifile)
	for surf, (pshort, levs, vars) in base_data.items():
		if os.path.exists(obfile % pshort):
			print('Reading cached data for level %s' % pshort)
			fc = np.load(obfile % pshort)
			bdat[pshort] = dict(fc)
			fc.close()
		
		else:
			bdat[pshort] = {}
			ckwargs = {'nameOfFirstFixedSurface': surf}
			if not levs:
				raise NotImplementedError('None should mean all available levels.')
			i = 0
			zlen = len(levs)
			for lev in levs:
				print('Reading level %s:%s' % (pshort,lev))
				ckwargs['scaledValueOfFirstFixedSurface'] = lev
				for var in vars:
					msg = f.select(shortNameECMF=var, **ckwargs)
					if len(msg) > 1:
						raise ValueError('Got more than one message matching the query!')
					msg = msg[0]

					if i == 0:
						bdat[pshort][var] = np.empty((1, zlen, ) + msg.values.shape)
					dat = np.array(msg.values.data)
					dat[dat == msg.missingValue] = np.nan
					bdat[pshort][var][0,i,:,:] = dat
				
				if i == 0:
					bdat[pshort]['lvl'] = np.empty((zlen,))
				bdat[pshort]['lvl'][i] = msg.scaledValueOfFirstFixedSurface
				
				i += 1
			
			size = 0
			recs = 0
			for var in bdat[pshort]:
				bdat[pshort][var] = bdat[pshort][var].astype('f4')
				size += bdat[pshort][var].size * 4
				recs += bdat[pshort][var].size/720/361
			print('Saving basic data on %s level (%d fields, %.2fM uncompressed)' % (pshort, recs, size/1024**2))
			np.savez_compressed(obfile % pshort, **bdat[pshort])

		ddat[pshort] = {}
		
	f.close()

	# Calculate reference diagnostics
	for pshort, diagnostics_ in diagnostics.items():
		calc(bdat[pshort], ddat[pshort], grid, diagnostics_, slicebylev=False)

	# Save results
	for pshort in ddat.keys():
		size = 0
		recs = 0
		for var in ddat[pshort]:
			ddat[pshort][var] = ddat[pshort][var].astype('f4')
			size += ddat[pshort][var].size * 4
			recs += ddat[pshort][var].size/720/361
		print('Saving derived data on %s level (%d fields, %.2fM uncompressed)' % (pshort, recs, size/1024**2) )
		np.savez_compressed(odfile % pshort, **ddat[pshort])
	
# C'est le fin
