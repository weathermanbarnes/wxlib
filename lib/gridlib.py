#!/usr/bin/env python
# -*- encoding: utf-8


''' Obtain meta information about data

The meta-information currently mainly concerns the grid information that is 
required for plotting.

This information is either extracted from a netCDF file object, or taken from
a "static.npz" file that contains the pertinent information for a given data set.
'''


import copy
import numpy as np
import netCDF4 as nc
from mpl_toolkits.basemap.pyproj import Proj 		# for rotpole projection
from datetime import datetime as dt, timedelta as td

from . import derivatives


class grid(object):
	''' Underlying object defining the grid API and basic grid types '''

	def __init__(self):
		self.gridtype = 'unset'
		self.x = None
		self.x_name = None
		self.y = None
		self.y_name = None
		self.z = None
		self.z_name = None
		self.t = None
		self.t_name = None
		self.nx = 0
		self.ny = 0
		self.nz = 0
		self.nt = 0

		self._init_grid()
		self.__build_grid()
		
		# Basemap won't work with Fortran-aligned arrays
		self.x = np.ascontiguousarray(self.x)
		self.y = np.ascontiguousarray(self.y)

		return
	
	# To be overwritten by derived classes
	def _init_grid(self):
		pass

	def _calc_dx_dy_latlon(self):
		self.dx = np.ones((self.ny, self.nx))*111111.111111
		self.dy = np.ones((self.ny, self.nx))*111111.111111

		dlon = np.ones(self.dx.shape)
		dlon[:,1:-1] = self.x[np.newaxis,2:]-self.x[np.newaxis,:-2]
		dlon[:,0] = self.x[np.newaxis,1]-self.x[np.newaxis,-1]
		if np.any(np.abs(dlon[:,1]/dlon[:,0]-1)  > 1.0e-3):
			self.cyclic_ew = False
			dlon[:,0] = 2.0*(self.x[1]-self.x[0])
			dlon[:,-1] = 2.0*(self.x[-1]-self.x[-2])
		else:
			dlon[:,-1] = self.x[np.newaxis,0]-self.x[np.newaxis,-2]
		dlon[dlon > 180] -= 360.0
		dlon[dlon < -180] += 360.0
		self.dx *= dlon * np.cos(np.pi/180.0*self.y[:,np.newaxis])

		self.dy[1:-1,:] *= self.y[2:,np.newaxis]-self.y[:-2,np.newaxis]
		self.dy[ 0,:] *= 2.0*(self.y[ 1]-self.y[ 0])
		self.dy[-1,:] *= 2.0*(self.y[-1]-self.y[-2])

		return

	# Building the grid on top of the established axes
	def __build_grid(self):
		# rather obsolete ...

		return

	def _parse_time(self):
		tusplit = self.t_unit.split()
		if len(tusplit) > 3 and tusplit[1] == 'since':
			facs = {'seconds': 1, 'minutes': 60, 'hours': 3600, 'days': 86400}
			if tusplit[0] not in facs:
				self.t_parsed = None
			else:
				formats = ['%Y-%m-%d %H:%M:0.0', '%Y-%m-%d %H:%M:00', '%Y-%m-%d %H:%M', ]
				self.t_interval_unit = facs[tusplit[0]]
				for fmt in formats:
					try:
						self.t_epoch = dt.strptime(' '.join(tusplit[2:4]), fmt)
					except ValueError:
						self.t_epoch = None

					if not type(self.t_epoch) == type(None):
						break

				if type(self.t_epoch) == type(None):
					self.t_epoch = dt(1,1,1,0,0)
					
				self.t_parsed = np.array([self.t_epoch + td(0, self.t_interval_unit*int(ts)) for ts in self.t])

		else:
			self.t_parsed = None
		
		return

	rebuild_grid = __build_grid

	def new_time(self, dates):
		''' Makes a copy of itself with a new time axis

		Parameters
		----------
		dates : list of datetime
		    Dates defining the new time axis

		Returns
		-------
		    Copy of the grid object with a new time axis
		'''
		
		if type(self.t_parsed) == type(None):
			raise TypeError('No (parsable) time axis to be replaced in this grid object!')

		cpy = copy.copy(self)
		cpy.t_parsed = dates
		cpy.t = np.array([ (date - cpy.t_epoch).total_seconds()/float(cpy.t_interval_unit) 
			for date in dates ])

		return cpy

	def unrotate_vector(self, u, v):
		''' Express vector components expressed in rotated grid in unrotated grid

		Only active if grid instance represents a rotated grid. Otherwise u and v a returned unchanged.

		Parameters
		----------
		u : np.ndarray with 2 or more dimensions
		    x-component of vector in rotated coordinates
		v : np.ndarray with 2 or more dimensions
		    x-component of vector in rotated coordinates

		Returns
		----------
		np.ndarray with 2 or more dimensions
		    x-component of vector in unrotated coordinates
		np.ndarray with 2 or more dimensions
		    x-component of vector in unrotated coordinates
		'''

		if not hasattr(self, 'local_Nx'):
			return u, v
		
		ur = u * self.local_Ny - v * self.local_Nx
		vr = u * self.local_Nx + v * self.local_Ny

		return ur, vr




# Construct the grid based on the grid information in a nc (netcdf) file
class grid_by_nc(grid):
	''' Extract the relevant information from a given netCDF file object '''

	X_NAMES = ['lon', 'longitude', 'west_east', 'west_east_stag', 'x', 'x_1']
	X_NAME_BEGINSWITH = ['rlon', 'dimx', ]
	Y_NAMES = ['lat', 'latitude', 'south_north', 'south_north_stag', 'y', 'y_1']
	Y_NAME_BEGINSWITH = ['rlat', 'dimy', ]
	Z_NAMES = ['level', 'bottom_top', 'bottom_top_stag', 'z', 'z_1', 'alt']
	Z_NAME_BEGINSWITH = ['lev', 'dimz', ]
	T_NAMES = ['time', 'Time']

	ROT_POLES = {
		'rotated_pole': ('grid_north_pole_longitude', 'grid_north_pole_latitude'),  # name convention in NORA10
	}

	def __init__(self, ncfile, ncvar=None):
		self.f = ncfile
		self.v = ncvar
		self.oro = None

		grid.__init__(self)
		
		return

	# Skims through the netcdf file looking for the type of the x and y axis
	def _init_grid(self):
		def matches(match, names, begins):
			for name in names:
				if name == match: 
					return True
			for begin in begins:
				if begin == match[:len(begin)]:
					return True
			
			return False
	
		# Part 1: Looking for suitable axis
		# 1. Using the given variable
		if self.v:
			for d in self.v.dimensions:
				if matches(d, self.X_NAMES, self.X_NAME_BEGINSWITH):
					if self.x_name:
						raise ValueError('Found several possible x-axes (using variable)')
					self.x_name = d
				if matches(d, self.Y_NAMES, self.Y_NAME_BEGINSWITH):
					if self.y:
						raise ValueError('Found several possible y-axes (using variable)')
					self.y_name = d
				if matches(d, self.Z_NAMES, self.Z_NAME_BEGINSWITH):
					if self.z_name:
						raise ValueError('Found several possible z-axes (using variable)')
					self.z_name = d
				if d in self.T_NAMES:
					if self.t_name:
						raise ValueError('Found several possible t-axes (using variable)')
					self.t_name = d

		if not self.x_name or not self.y_name or not self.z_name or not self.t_name:
			self.x_name = None
			self.y_name = None
			self.z_name = None
			self.t_name = None

			for d in self.f.dimensions:
				if matches(d, self.X_NAMES, self.X_NAME_BEGINSWITH):
					if self.x_name:
						raise ValueError('Found several possible x-axes (using file)')
					self.x_name = d
				if matches(d, self.Y_NAMES, self.Y_NAME_BEGINSWITH):
					if self.y:
						raise ValueError('Found several possible y-axes (using file)')
					self.y_name = d
				if matches(d, self.Z_NAMES, self.Z_NAME_BEGINSWITH):
					if self.z_name:
						raise ValueError('Found several possible z-axes (using file)')
					self.z_name = d
				if d in self.T_NAMES:
					if self.t_name:
						raise ValueError('Found several possible t-axes (using file)')
					self.t_name = d
		
		if not self.x_name:
			raise ValueError('No x-axis found')
		if not self.y_name:
			raise ValueError('No y-axis found')

		# Part 2: Determining type of axis
		self.gridtype = None
		self.cyclic_ew = False
		self.cyclic_ns = False
		
		try:
			self.x_unit = self.f.variables[self.x_name].units
		except (KeyError, AttributeError):
			self.x_unit = '1'
		try: 
			self.y_unit = self.f.variables[self.y_name].units
		except (KeyError, AttributeError): 
			self.y_unit = '1'
		if self.z_name:
			try:
				self.z_unit = self.f.variables[self.z_name].units
			except (KeyError, AttributeError):
				self.z_unit = '1'
		else:
			self.z_unit = None
		if self.t_name:
			try:
				self.t_unit = self.f.variables[self.t_name].units
			except (KeyError, AttributeError):
				self.t_unit = '1'
		else:
			self.t_unit = None
		
		self.nx = len(self.f.dimensions[self.x_name])
		self.ny = len(self.f.dimensions[self.y_name])

		if self.x_unit == 'degrees_E' and self.y_unit == 'degrees_N':
			self.gridtype = 'latlon'
			self.cyclic_ew = True
			self.x = self.f.variables[self.x_name][::]
			self.y = self.f.variables[self.y_name][::]
			self.x_name = 'longitude'
			self.y_name = 'latitude'
		elif self.x_unit == 'degrees_east' and self.y_unit == 'degrees_north':
			self.gridtype = 'latlon'
			self.cyclic_ew = True
			self.x = self.f.variables[self.x_name][::]
			self.y = self.f.variables[self.y_name][::]
			self.x_name = 'longitude'
			self.y_name = 'latitude'
		elif self.x_unit in ['degrees', 'degree'] and self.y_unit in ['degrees', 'degree']:
			self.gridtype = 'latlon'
			self.cyclic_ew = True
			self.x = self.f.variables[self.x_name][::]
			self.y = self.f.variables[self.y_name][::]
			self.x_name = 'longitude'
			self.y_name = 'latitude'
		elif self.x_unit == '1' and self.y_unit == '1':
			if 'OUTPUT FROM WRF' in getattr(self.f, 'TITLE', ''):
				if self.f.GRIDTYPE == 'C':
					self.gridtype = 'cartesian'
					self.x = np.arange(self.nx)*self.f.DX
					self.y = np.arange(self.ny)*self.f.DY
					# Just assuming that WRF is being sensible.
					self.x_unit = 'm'
					self.y_unit = 'm'
					self.x_name = 'x'
					self.y_name = 'y'

				else:
					raise NotImplementedError('Unknown WRF gridtype `%s\'' % self.f.GRIDTYPE)
			else:
				self.gridtype = 'idx'
				self.x = np.arange(self.nx)
				self.y = np.arange(self.ny)
				self.x_name = 'xidx'
				self.y_name = 'yidx'
		elif self.x_unit == 'm' and self.y_unit == 'm':
			self.gridtype = 'cartesian'
			self.x = self.f.variables[self.x_name][::]
			self.y = self.f.variables[self.y_name][::]
			self.x_name = 'x'
			self.y_name = 'y'

		else:
			raise NotImplementedError('(Yet) Unknown grid type with units (%s/%s)' % (self.x_unit, self.y_unit))
			
		if self.gridtype == 'latlon':
			self._calc_dx_dy_latlon()
			self.x = np.tile(self.x, (self.ny,1))
			self.y = np.tile(self.y, (self.nx,1)).T

		elif self.gridtype == 'idx':
			self.dx = np.ones((self.ny, self.nx))*2
			self.x = np.tile(self.x, (self.ny,1))
			self.y = np.tile(self.y, (self.nx,1)).T
			self.dy = np.ones((self.ny, self.nx))*2

		elif self.gridtype == 'cartesian':
			self.x = np.tile(self.x, (self.ny,1))
			self.y = np.tile(self.y, (self.nx,1)).T
			self.dx = np.empty(self.x.shape)
			self.dy = np.empty(self.y.shape)
			self.dx[:,1:-1] = self.x[:,2:] - self.x[:,:-2]
			self.dy[1:-1,:] = self.y[2:,:] - self.y[:-2,:]

		else:
			raise NotImplementedError('(Yet) Unknown grid type "%s"' % self.gridtype)

		self.rotated = False
		if self.gridtype == 'latlon':
			for var in self.f.variables:
				if var in self.ROT_POLES:
					rot_nplon_name, rot_nplat_name = self.ROT_POLES[var]
					rot_nplon = getattr(self.f.variables[var], rot_nplon_name) 
					rot_nplat = getattr(self.f.variables[var], rot_nplat_name)
					m = Proj(proj='ob_tran', o_proj='latlon', 
							o_lon_p=rot_nplon,
							o_lat_p=rot_nplat, 
							lon_0=180,
					)
					self.x, self.y = m(np.ascontiguousarray(self.x), 
							np.ascontiguousarray(self.y) )
					self.x *= 180.0/np.pi
					self.y *= 180.0/np.pi
					self.rotated = True
					self.local_Nx, self.local_Ny = derivatives.grad(self.y[np.newaxis,:,:], self.dx, self.dy)
					self.local_Nx, self.local_Ny = self.local_Nx.squeeze(), self.local_Ny.squeeze()
					absgrad = np.sqrt(self.local_Nx**2 + self.local_Ny**2)
					self.local_Nx /= absgrad
					self.local_Ny /= absgrad
					break # don't rotate more than once!

		if self.z_name:
			self.nz = len(self.f.dimensions[self.z_name])
			if self.z_name in self.f.variables:
				self.z  = self.f.variables[self.z_name][::]
			else:
				self.z = np.arange(self.nz)
		if self.t_name:
			self.nt = len(self.f.dimensions[self.t_name])
			if not self.nt and not self.v: 
				raise RuntimeError('grid_by_nc needs one specific variable for extracing the length of the netcdf-unlimited time dimension')
			elif not self.nt:
				timedim = self.v.dimensions.index(self.t_name)
				self.nt = self.v.shape[timedim]
			if self.t_name in self.f.variables:
				self.t = self.f.variables[self.t_name][::]
			else:
				self.t = np.arange(self.nt)
			
			# Try to parse the time axis into datetime objects
			t = self.f.variables[self.t_name] 
			if hasattr(t, 'units'):
				self.t_parsed = nc.num2date(t[:], units=t.units, calendar=getattr(t, 'calendar', 'standard'))

		return


# Construct the grid based on the information in a static.npz file
class grid_by_static(grid):
	''' Read relevant information from a static file, usually called "static.npz" '''

	def __init__(self, staticfile):
		self.f = staticfile

		grid.__init__(self)

		return
	

	def _init_grid(self):
		if 'lat' in self.f.files and 'lon' in self.f.files:
			self.gridtype = 'latlon'
			self.cyclic_ew = True
			self.oro = None
			
			self.nx = self.f['lon'].shape[0]
			self.ny = self.f['lat'].shape[0]
			self.x  = self.f['lon'][:]
			self.y  = self.f['lat'][:]
			self.x_name = 'longitude'
			self.y_name = 'latitude'
			self.x_unit = 'degrees_east'
			self.y_unit = 'degrees_north'

			self._calc_dx_dy_latlon()
		else:
			raise NotImplementedError('(Yet) Unknown grid type using the variables ' % str(self.f.files))

		self.x = np.tile(self.x, (self.ny,1))
		self.y = np.tile(self.y, (self.nx,1)).T

		return


#
