#!/usr/bin/env python
# -*- encoding: utf-8


''' Obtain meta information about data

The meta-information currently mainly concerns the grid information that is 
required for plotting.

This information is either extracted from a netCDF file object, or taken from
a "static.npz" file that contains the pertinent information for a given data set.
'''


import warnings
import copy
import numpy as np
import netCDF4 as nc

# Required only for rotpole projection, so this should not be a hard dependency
NO_PROJ_AVAILABLE = False
try:
    from mpl_toolkits.basemap.pyproj import Proj
except:
    try: 
        from pyproj import Proj
    except:
        NO_PROJ_AVAILABLE = True

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

        if len(self.x.shape) == 1:
            slc = np.newaxis
        else:
            slc = slice(None)

        dlon = np.ones(self.dx.shape)
        dlon[:,1:-1] = (self.x[slc,2:]-self.x[slc,:-2]) % 360.0
        dlon[:,0] = (self.x[slc,1]-self.x[slc,-1]) % 360.0
        if np.any(np.abs(dlon[:,1]/dlon[:,0]-1)  > 1.0e-3):
            self.cyclic_ew = False
            dlon[:,0] = 2.0*(self.x[slc,1]-self.x[slc,0])
            dlon[:,-1] = 2.0*(self.x[slc,-1]-self.x[slc,-2])
        else:
            dlon[:,-1] = (self.x[slc,0]-self.x[slc,-2]) % 360.0
        dlon[dlon > 180] -= 360.0
        dlon[dlon < -180] += 360.0
        self.dx *= dlon * np.cos(np.pi/180.0*self.y[:,slc])

        self.dy[1:-1,:] *= self.y[2:,slc]-self.y[:-2,slc]
        self.dy[ 0,:] *= 2.0*(self.y[ 1,slc]-self.y[ 0,slc])
        self.dy[-1,:] *= 2.0*(self.y[-1,slc]-self.y[-2,slc])

        return
    
    def _calc_dx_dy_cartesian(self):
        self.dx = np.empty(self.x.shape)
        self.dx[:,1:-1] = self.x[:,2:] - self.x[:,:-2]
        self.dx[:,0] = (self.x[:,1] - self.x[:,0])*2
        self.dx[:,-1] = (self.x[:,-2] - self.x[:,-1])*2

        self.dy = np.empty(self.y.shape)
        self.dy[1:-1,:] = self.y[2:,:] - self.y[:-2,:]
        self.dy[0,:] = (self.y[1,:] - self.y[0,:])*2
        self.dy[-1,:] = (self.y[-1,:] - self.y[-2,:])*2

        return

    # Building the grid on top of the established axes
    def __build_grid(self):
        # rather obsolete ...

        return

    def _parse_time_unit(self):
        tusplit = self.t_unit.split()
        if len(tusplit) > 3 and tusplit[1] == 'since':
            facs = {'seconds': 1, 'minutes': 60, 'hours': 3600, 'days': 86400}
            if tusplit[0] in facs:
                formats = ['%Y-%m-%d %H:%M:00.0', '%Y-%m-%d %H:%M:0.0', '%Y-%m-%d %H:%M:00', '%Y-%m-%d %H:%M', ]
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
        
        if not hasattr(self, 't_epoch') or not hasattr(self, 't_interval_unit'):
            raise TypeError('Need a reference date and time interval!')

        if type(self.t_name) == type(None):
            self.t_name = 'time'

        cpy = copy.copy(self)
        cpy.t_parsed = dates
        cpy.t = np.array([ (date - cpy.t_epoch).total_seconds()/float(cpy.t_interval_unit) 
            for date in dates ])

        return cpy

    def new_plev(self, plev):
        cpy = copy.copy(self)

        if type(plev) == type(None):
            cpy.z_name = ''
            cpy.z_unit = ''
            cpy.z = []
        elif plev == 'sfc':
            cpy.z_name = None
            cpy.z = None
            cpy.z_unit = None
        elif plev[:2] == 'pv':
            cpy.z_name = 'lev'
            cpy.z_unit = 'PVU'
            cpy.z = [int(plev[2:])/1000, ]
        elif plev[:2] == 'pt':
            cpy.z_name = 'lev'
            cpy.z_unit = 'K'
            cpy.z = [int(plev[2:]), ]
        else:
            cpy.z_name = 'lev'
            cpy.z_unit = 'Pa'
            cpy.z = [int(plev)*100, ]

        return cpy

    def prepare_rotation(self):
        self.rotated = True
        self.local_Nx, self.local_Ny = derivatives.grad(self.y[np.newaxis,:,:], self.dx, self.dy)
        self.local_Nx, self.local_Ny = self.local_Nx.squeeze(), self.local_Ny.squeeze()
        absgrad = np.sqrt(self.local_Nx**2 + self.local_Ny**2)
        self.local_Nx /= absgrad
        self.local_Ny /= absgrad

        return


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

    X_NAMES = ['lon', 'longitude', 'west_east', 'west_east_stag', 'x', 'x_1', 'x (stag)']
    X_NAME_BEGINSWITH = ['rlon', 'srlon', 'dimx', ]
    Y_NAMES = ['lat', 'latitude', 'south_north', 'south_north_stag', 'y', 'y_1', 'y (stag)']
    Y_NAME_BEGINSWITH = ['rlat', 'srlat', 'dimy', ]
    Z_NAMES = ['level', 'bottom_top', 'bottom_top_stag', 'z', 'z_1', 'alt', 'z (stag)']
    Z_NAME_BEGINSWITH = ['lev', 'dimz', ]
    T_NAMES = ['time', 'Time']

    X_IRREGULAR_NAMES = X_NAMES
    X_IRREGULAR_NAME_BEGINSWITH = X_NAME_BEGINSWITH
    Y_IRREGULAR_NAMES = Y_NAMES
    Y_IRREGULAR_NAME_BEGINSWITH = Y_NAME_BEGINSWITH

    ROT_POLES = {
        'rotated_pole': ('grid_north_pole_longitude', 'grid_north_pole_latitude'),  # name convention in NORA10, COSMO and dynlib
        'projection_3': ('grid_north_pole_longitude', 'grid_north_pole_latitude'),  # name in NORA10 altitude file
    }
    VAR_PROJ = ['projection_lambert', ] # Marks the grid projection for AROME

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
                if name == match.lower(): 
                    return True
            for begin in begins:
                if begin == match.lower()[:len(begin)]:
                    return True
            
            return False
    
        # Part 1: Looking for suitable axis
        # 1. Using the given variable
        if self.v:
            for d in self.v.dimensions:
                if matches(d, self.X_NAMES, self.X_NAME_BEGINSWITH):
                    if self.x_name:
                        raise ValueError('Found several possible x-axes: %s, %s (using variable)' % (self.x_name, d))
                    self.x_name = d
                if matches(d, self.Y_NAMES, self.Y_NAME_BEGINSWITH):
                    if self.y:
                        raise ValueError('Found several possible y-axes: %s, %s (using variable)' % (self.y_name, d))
                    self.y_name = d
                if matches(d, self.Z_NAMES, self.Z_NAME_BEGINSWITH):
                    if self.z_name:
                        raise ValueError('Found several possible z-axes: %s, %s (using variable)' % (self.z_name, d))
                    self.z_name = d
                if d in self.T_NAMES:
                    if self.t_name:
                        raise ValueError('Found several possible t-axes: %s, %s (using variable)' % (self.t_name, d))
                    self.t_name = d
            
            # In case of irregular grids, coordinates might be given as separate coordinate variables
            if not self.x_name or not self.y_name:
                for cv in self.f.variables:
                    if matches(cv, self.X_IRREGULAR_NAMES, self.X_IRREGULAR_NAME_BEGINSWITH):
                        if cv == self.x_name:
                            continue
                        if self.x_name:
                            raise ValueError('Found several possible x-axes: %s, %s (scanning for separate coordinate variables)' % (self.x_name, cv))
                        if self.f.variables[cv].dimensions == self.v.dimensions[-2:]:
                            self.x_name = cv
                    if matches(cv, self.Y_IRREGULAR_NAMES, self.Y_IRREGULAR_NAME_BEGINSWITH):
                        if cv == self.y_name:
                            continue
                        if self.y_name:
                            raise ValueError('Found several possible y-axes: %s, %s (scanning for separate coordinate variables)' % (self.y_name, cv))
                        if self.f.variables[cv].dimensions == self.v.dimensions[-2:]:
                            self.y_name = cv

        if not self.x_name or not self.y_name:
            self.x_name = None
            self.y_name = None
            self.z_name = None
            self.t_name = None

            give_up_z = False
            give_up_t = False

            for d in self.f.dimensions:
                if matches(d, self.X_NAMES, self.X_NAME_BEGINSWITH):
                    if self.x_name:
                        raise ValueError('Found several possible x-axes: %s, %s (using file)' % (self.x_name, d))
                    self.x_name = d
                if matches(d, self.Y_NAMES, self.Y_NAME_BEGINSWITH):
                    if self.y:
                        raise ValueError('Found several possible y-axes: %s, %s (using file)' % (self.y_name, d))
                    self.y_name = d
                if matches(d, self.Z_NAMES, self.Z_NAME_BEGINSWITH) and not give_up_z:
                    if self.z_name:
                        warnings.warn('Found several possible z-axes: %s, %s (using file)' % (self.z_name, d))
                        self.z_name = None
                        give_up_z = True
                    self.z_name = d
                if d in self.T_NAMES and not give_up_t:
                    if self.t_name:
                        warnings.warn('Found several possible t-axes: %s, %s (using file)' % (self.t_name, d))
                        self.t_name = None
                        give_up_t = True
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
        
        self.nx = self.f.dimensions[self.x_name].size
        self.ny = self.f.dimensions[self.y_name].size

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
        elif self.x_unit == 'degrees east' and self.y_unit == 'degrees north':
            self.gridtype = 'latlon'
            self.cyclic_ew = True
            self.x = self.f.variables[self.x_name][::]
            self.y = self.f.variables[self.y_name][::]
            self.x_name = 'longitude'
            self.y_name = 'latitude'
        elif self.x_unit == 'degree_east' and self.y_unit == 'degree_north':
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
        elif self.x_unit in ['m', 'km'] and self.y_unit in ['km', 'm']:
            self.gridtype = 'cartesian'
            self.x = self.f.variables[self.x_name][::]
            self.y = self.f.variables[self.y_name][::]
            self.x_name = 'x'
            self.y_name = 'y'

        else:
            raise NotImplementedError('(Yet) Unknown grid type with units (%s/%s)' % (self.x_unit, self.y_unit))
            
        if self.gridtype == 'latlon':
            self._calc_dx_dy_latlon()
            if len(self.x.shape) == 1:
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
            self._calc_dx_dy_cartesian()

        else:
            raise NotImplementedError('(Yet) Unknown grid type "%s"' % self.gridtype)

        self.rotated = False
        if self.gridtype == 'latlon':
            for var in self.f.variables:
                if var in self.ROT_POLES:
                    if NO_PROJ_AVAILABLE:
                        raise ValueError('For handling rotated-pole grids, a Proj module must be available')
                    
                    rot_nplon_name, rot_nplat_name = self.ROT_POLES[var]
                    rot_nplon = getattr(self.f.variables[var], rot_nplon_name) 
                    rot_nplat = getattr(self.f.variables[var], rot_nplat_name)
                    m = Proj(proj='ob_tran', o_proj='latlon', 
                            o_lon_p=rot_nplon,
                            o_lat_p=rot_nplat, 
                            lon_0=180,
                    )
                    self.rot_np = (rot_nplat, rot_nplon)
                    self.rot_x_name, self.rot_y_name = 'rlon', 'rlat'
                    self.rot_x_longname, self.rot_y_longname = 'rotated_longitude', 'rotated_latitude'
                    self.rot_x, self.rot_y = self.x, self.y
                    self.x, self.y = m(np.ascontiguousarray(self.x), 
                            np.ascontiguousarray(self.y) )
                    self.x *= 180.0/np.pi
                    self.y *= 180.0/np.pi
                    self.prepare_rotation()
                    break # don't rotate more than once!
                if var in self.VAR_PROJ:
                    self.prepare_rotation()
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
                self.nt = 0
            elif not self.nt:
                if self.t_name in self.v.dimensions:
                    timedim = self.v.dimensions.index(self.t_name)
                    self.nt = self.v.shape[timedim]
                else:
                    self.nt = 0
            if self.t_name in self.f.variables:
                self.t = self.f.variables[self.t_name][::]

                # Try to parse the time axis into datetime objects
                t = self.f.variables[self.t_name] 
                if hasattr(t, 'units'):
                    self.t_unit = t.units
                    self._parse_time_unit() # set self.t_epoch and self.t_interval_unit
                    
                    # A unit unknown by cftime as of cftime 1.5.1.1, encountered in OpenIFS output
                    if t.units == 'day as %Y%m%d.%f':
                        tvals = t[:]
                        self.t_parsed = [dt(int(tval/10000), int((tval / 100) % 100), int(tval % 100)) + td(tval % 1) 
                                for tval in tvals]

                    else:
                        self.t_parsed = nc.num2date(t[:], units=t.units, calendar=getattr(t, 'calendar', 'standard'))
            else:
                self.t = np.arange(self.nt)
                self.t_parsed = None

        return





# Construct the grid based on given lat/lon arrays
class grid_by_latlon(grid):
    ''' Build the relevant information from given lats and lons '''

    def __init__(self, lats, lons):
        self.lon = lons
        self.lat = lats
        
        grid.__init__(self)

        return

    # Skims through the netcdf file looking for the type of the x and y axis
    def _init_grid(self):
        self.gridtype = 'latlon'
        self.cyclic_ew = True
        self.cyclic_ns = False
        
        self.x = self.lon
        self.y = self.lat
        self.x_name = 'longitude'
        self.y_name = 'latitude'
        self.x_unit = 'degrees_east'
        self.y_unit = 'degrees_north'

        self.ny, self.nx = self.lon.shape
        
        self._calc_dx_dy_latlon()
        if len(self.x.shape) == 1:
            self.x = np.tile(self.x, (self.ny,1))
            self.y = np.tile(self.y, (self.nx,1)).T

        self.rotated = False

        return


# Construct the grid based on given x / y arrays
class grid_by_xy(grid):
    ''' Build the relevant information from given lats and lons '''

    def __init__(self, x, y, dx=None, dy=None, lats=None, lons=None, oro=None):
        self.__init_args = x, y, dx, dy, lats, lons
        grid.__init__(self)

        if not type(oro) == type(None):
            self.oro = oro

        return

    # Skims through the netcdf file looking for the type of the x and y axis
    def _init_grid(self):
        self.gridtype = 'cartesian'
        self.cyclic_ew = False
        self.cyclic_ns = False
        
        x, y, dx, dy, lats, lons = self.__init_args

        self.x = x
        self.y = y

        if not type(dx) == type(None):
            self.dx = dx
            self.dy = dy
        else:
            self._calc_dx_dy_cartesian()

        if not type(lats) == type(None):
            self.lon = lons
            self.lat = lats
        
        self.x_name = 'x'
        self.y_name = 'y'
        self.x_unit = 'm'
        self.y_unit = 'm'

        self.ny, self.nx = self.x.shape
        
        self.rotated = False

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
