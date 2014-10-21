#!/usr/bin/env python
# -*- encoding: utf-8
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  		DynLib -- Grid definitions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
import math
import numpy as np
from datetime import datetime as dt, timedelta as td


class grid(object):
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

		return
	
	# To be overwritten by derived classes
	def _init_grid(self):
		pass

	# Building the grid on top of the established axes
	def __build_grid(self):
		if self.gridtype == 'latlon':
			self.dx = np.ones((self.ny, self.nx))*111111.111111
			self.dy = np.ones((self.ny, self.nx))*111111.111111
			for xidx in range(self.nx):
				for yidx in range(1,self.ny-1):
					dlon = self.x[(xidx+1)%self.nx]-self.x[(xidx-1)%self.nx]
					if dlon > 180:
						dlon -= 360
					elif dlon < -180:
						dlon += 360
					self.dx[yidx,xidx] *= dlon*math.cos(math.pi/180.0*self.y[yidx])
			for yidx in range(1,self.ny-1):
				self.dy[yidx,:] *= self.y[yidx+1]-self.y[yidx-1]
			self.dy[ 0,:] *= 2.0*(self.y[ 1]-self.y[ 0])
			self.dy[-1,:] *= 2.0*(self.y[-1]-self.y[-2])

			self.x = np.tile(self.x, (self.ny,1))
			self.y = np.tile(self.y, (self.nx,1)).T

		elif self.gridtype == 'idx':
			self.dx = np.ones((self.ny, self.nx))*2
			self.dy = np.ones((self.ny, self.nx))*2

			self.x = np.tile(self.x, (self.ny,1))
			self.y = np.tile(self.y, (self.nx,1)).T
		elif self.gridtype == 'cartesian':
			self.x = np.tile(self.x, (self.ny,1))
			self.y = np.tile(self.y, (self.nx,1)).T

			self.dx = np.empty(self.x.shape)
			self.dy = np.empty(self.y.shape)
			self.dx[:,1:-1] = self.x[:,2:] - self.x[:,:-2]
			self.dy[1:-1,:] = self.y[2:,:] - self.y[:-2,:]

		else:
			raise NotImplementedError, '(Yet) Unknown grid type "%s"' % self.gridtype

		if type(self.x) == np.ndarray: self.x = self.x[self.cut]
		if type(self.y) == np.ndarray: self.y = self.y[self.cut]
		if type(self.z) == np.ndarray: self.z = self.z[self.cut]

		if type(self.dx) == np.ndarray: self.dx = self.dx[self.cut]
		if type(self.dy) == np.ndarray: self.dy = self.dy[self.cut]

		if type(self.oro) == np.ndarray: self.oro = self.oro[self.cut]

		return


# Construct the grid based on the grid information in a nc (netcdf) file
class grid_by_nc(grid):
	X_NAMES = ['lon', 'longitude', 'west_east', 'west_east_stag', 'x']
	Y_NAMES = ['lat', 'latitude', 'south_north', 'south_north_stag', 'y']
	Z_NAMES = ['level', 'bottom_top', 'bottom_top_stag', 'z']
	T_NAMES = ['time', 'Time']

	def __init__(self, ncfile, ncvar=None, cut=slice(None)):
		self.f = ncfile
		self.v = ncvar
		self.cut = cut
		self.oro = None

		grid.__init__(self)
		
		return

	# Skims through the netcdf file looking for the type of the x and y axis
	def _init_grid(self):
		# Part 1: Looking for suitable axis
		if not self.v:
			lookout = self.f
		else:
			lookout = self.v
			
		for d in lookout.dimensions:
			if d in self.X_NAMES:
				if self.x_name:
					raise ValueError, 'Found several possible x-axes (did you set a variable?)'
				self.x_name = d
			if d in self.Y_NAMES:
				if self.y:
					raise ValueError, 'Found several possible y-axes (did you set a variable?)'
				self.y_name = d
			if d in self.Z_NAMES:
				if self.z_name:
					raise ValueError, 'Found several possible z-axes (did you set a variable?)'
				self.z_name = d
			if d in self.T_NAMES:
				if self.t_name:
					raise ValueError, 'Found several possible t-axes (did you set a variable?)'
				self.t_name = d
		
		if not self.x_name:
			raise ValueError, 'No x-axis found'
		if not self.y_name:
			raise ValueError, 'No y-axis found'

		# Part 2: Determining type of axis
		self.gridtype = None
		self.cyclic_ew = False
		self.cyclic_ns = False
		
		try:
			self.x_unit = self.f.variables[self.x_name].units
		except KeyError:
			self.x_unit = '1'
		try: 
			self.y_unit = self.f.variables[self.y_name].units
		except KeyError: 
			self.y_unit = '1'
		if self.z_name:
			try:
				self.z_unit = self.f.variables[self.z_name].units
			except KeyError:
				self.z_unit = '1'
		else:
			self.z_unit = None
		if self.t_name:
			try:
				self.t_unit = self.f.variables[self.t_name].units
			except KeyError:
				self.t_unit = '1'
		else:
			self.t_unit = None

		self.nx = self.f.dimensions[self.x_name]
		self.ny = self.f.dimensions[self.y_name]

		if self.x_unit == 'degrees_E' and self.y_unit == 'degrees_N':
			self.gridtype = 'latlon'
			self.cyclic_ew = True
			self.x = self.f.variables[self.x_name][::]
			self.y = self.f.variables[self.y_name][::]
		elif self.x_unit == 'degrees_east' and self.y_unit == 'degrees_north':
			self.gridtype = 'latlon'
			self.cyclic_ew = True
			self.x = self.f.variables[self.x_name][::]
			self.y = self.f.variables[self.y_name][::]
		elif self.x_unit == '1' and self.y_unit == '1':
			if 'TITLE' in self.f._attributes and 'OUTPUT FROM WRF' in self.f._attributes['TITLE']:
				if self.f._attributes['GRIDTYPE'] == 'C':
					self.gridtype = 'cartesian'
					self.x = np.arange(self.nx)*self.f._attributes['DX']
					self.y = np.arange(self.ny)*self.f._attributes['DY']
					# Just assuming that WRF is being sensible.
					self.x_unit = 'm'
					self.y_unit = 'm'

				else:
					raise NotImplementedError, 'Unknown WRF gridtype `%s\'' % self._attributes['GRIDTYPE']
			else:
				self.gridtype = 'idx'
				self.x = np.arange(self.nx)
				self.y = np.arange(self.ny)
		elif self.x_unit == 'm' and self.y_unit == 'm':
			self.gridtype = 'cartesian'
			self.x = self.f.variables[self.x_name][::]
			self.y = self.f.variables[self.y_name][::]

		else:
			raise NotImplementedError, '(Yet) Unknown grid type with units (%s/%s)' % (self.x_unit, self.y_unit)

		if self.z_name:
			self.nz = self.f.dimensions[self.z_name]
			if self.z_name in self.f.variables:
				self.z  = self.f.variables[self.z_name][::]
			else:
				self.z = np.arange(self.nz)
		if self.t_name:
			self.nt = self.f.dimensions[self.t_name]
			if not self.nt and not self.v: 
				raise RuntimeError, 'grid_by_nc needs one specific variable for extracing the length of the netcdf-unlimited time dimension'
			elif not self.nt:
				timedim = self.v.dimensions.index(self.t_name)
				self.nt = self.v.shape[timedim]
			if self.t_name in self.f.variables:
				self.t = self.f.variables[self.t_name][::]
			else:
				self.t = np.arange(self.nt)
			
			# Try to parse the time axis into datetime objects
			tusplit = self.t_unit.split()
			if tusplit[1] == 'since':
				facs = {'seconds': 1, 'minutes': 60, 'hours': 3600, 'days': 86400}
				if tusplit[0] not in facs:
					self.t_parsed = None
				else:
					formats = ['%Y-%m-%d %H:%M:0.0', '%Y-%m-%d %H:%M:00', ]
					fac = facs[tusplit[0]]
					for fmt in formats:
						try:
							start_dt = dt.strptime(' '.join(tusplit[2:4]), fmt)
						except ValueError:
							start_dt = None

						if not type(start_dt) == type(None):
							break

					if not type(start_dt) == type(None):
						start_dt = dt(1,1,1,0,0)
						
					self.t_parsed = [start_dt + td(0, fac*int(ts)) for ts in self.t]

			else:
				self.t_parsed = None

		return


# Construct the grid based on the information in a static.npz file
class grid_by_static(grid):
	def __init__(self, staticfile, cut=slice(None)):
		self.f = staticfile
		self.cut = cut

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
		else:
			raise NotImplementedError, '(Yet) Unknown grid type using the variables ' % str(self.f.files)

		return


#
