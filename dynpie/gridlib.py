#!/usr/bin/env python
# -*- encoding: utf-8
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  		DynLib -- Grid definitions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
import math
import numpy as np


class grid(object):
	def __init__(self):
		self.gridtype = 'unset'
		self.x = None
		self.y = None
		self.z = None
		self.t = None
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

		else:
			raise NotImplementedError, '(Yet) Unknown grid type "%s"' % self.gridtype

		return


# Construct the grid based on the grid information in a nc (netcdf) file
class grid_by_nc(grid):
	def __init__(self, ncfile, ncvar=None):
		self.f = ncfile
		self.v = ncvar

		grid.__init__(self)
		
		return

	# Skims through the netcdf file looking for the type of the x and y axis
	def _init_grid(self):
		# Part 1: Looking for suitable axis
		if not self.v:
			for d in self.f.dimensions:
				if d in ['lon','longitude']:
					if self.x:
						raise ValueError, 'Found several possible x-axes and no variable was given.'
					self.x = d
				if d in ['lat','latitude']:
					if self.y:
						raise ValueError, 'Found several possible y-axes and no variable was given.'
					self.y = d
				if d in ['level',]:
					if self.z:
						raise ValueError, 'Found several possible z-axes and no variable was given.'
					self.z = d
				if d in ['time',]:
					if self.t:
						raise ValueError, 'Found several possible z-axes and no variable was given.'
					self.t = d
		else:
			raise NotImplementedError
		
		if not self.x:
			raise ValueError, 'No x-axis found'
		if not self.y:
			raise ValueError, 'No y-axis found'

		# Part 2: Determining type of axis
		self.gridtype = None
		self.cyclic_ew = False
		self.cyclic_ns = False

		self.x_unit = self.f.variables[self.x].units
		self.y_unit = self.f.variables[self.y].units
		if self.z:
			self.z_unit = self.f.variables[self.z].units
		else:
			self.z_unit = None
		if self.t:
			self.t_unit = self.f.variables[self.t].units
		else:
			self.t_unit = None

		if self.x_unit == 'degrees_E' and self.y_unit == 'degrees_N':
			self.gridtype = 'latlon'
			self.cyclic_ew = True
		elif self.x_unit == 'degrees_east' and self.y_unit == 'degrees_north':
			self.gridtype = 'latlon'
			self.cyclic_ew = True
		else:
			raise NotImplementedError, '(Yet) Unknown grid type with units (%s/%s)' % (self.x_unit, self.y_unit)

		self.nx = self.f.variables[self.x].shape[0]
		self.ny = self.f.variables[self.y].shape[0]
		self.x  = self.f.variables[self.y][::]
		self.y  = self.f.variables[self.x][::]
		if self.z:
			self.nz = self.f.variables[self.z].shape[0]
			self.z  = self.f.variables[self.z][::]
		if self.t:
			self.nt = self.f.variables[self.t].shape[0]
			self.t  = self.f.variables[self.t][::]

		return


# Construct the grid based on the information in a static.npz file
class grid_by_static(grid):
	def __init__(self, staticfile):
		self.f = staticfile

		grid.__init__(self)

		return
	

	def _init_grid(self):
		if 'lat' in self.f.files and 'lon' in self.f.files:
			self.gridtype = 'latlon'
			self.cyclic_ew = True
			
			self.nx = self.f['lon'].shape[0]
			self.ny = self.f['lat'].shape[0]
			self.x  = self.f['lon'][:]
			self.y  = self.f['lat'][:]
		else:
			raise NotImplementedError, '(Yet) Unknown grid type using the variables ' % str(self.f.files)

		return


#
