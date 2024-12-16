##### Import Python Modules #####
import numpy as np
import pandas as pd

##### Import wxlib Modules #####
from utility import haversine
    
def grid_from_xarray(data,
                      lonname='longitude',latname='latitude',timename='time', **kwargs):

    ''' Create a grid dictionary from an xarray DataArray
   
    Parameters
    ----------
   
    data : xr.DataArray with shape (nz,ny,nx) and dtype float64

    Optional
    ----------
    lonname : name of the longitude variable in the xarray (default: longitude)
    latname : name of the latitude variable in the xarray (default: latitude)
    timename : name of the time variable in the xarray (default: time)
   
    Returns
    -------
    grid : dictionary of grid definitions required by detection

    grid dictionary includes:
        - longitude : 1D np.array of longitude
        - latitude : 1D np.array of latitude
        - time : 1D np.array of time (pd.datetime)
        - longitude2D : 2D np.array of longitude
        - latitude2D : 2D np.array of latitude
        - dx : x distance (m)
        - dy : y distance (m)
        - area : area of grid cell (m**2)
        - grid_distance : average grid distance over the domain (m)
    '''
    grid={}
    grid['longitude']=data[lonname].values
    grid['latitude']=data[latname].values
    grid['longitude2D'],grid['latitude2D']=np.meshgrid(data[lonname],data[latname])
    grid['time']=pd.to_datetime(data[timename])
    
    grid['dx'],grid['dy'],grid['area'],grid['grid_distance']=calc_grid_distance_area(grid['longitude2D'],grid['latitude2D'])

    return grid

def calc_grid_distance_area(lon,lat):
    """ Function to calculate grid parameters
        It uses haversine function to approximate distances
        It approximates the first row and column to the sencond
        because coordinates of grid cell center are assumed
        lat, lon: input coordinates(degrees) 2D [y,x] dimensions
        dx: distance (m)
        dy: distance (m)
        area: area of grid cell (m2)
        grid_distance: average grid distance over the domain (m)
    """
    dy = np.zeros(lon.shape)
    dx = np.zeros(lat.shape)

    dx[:,1:]=haversine(lon[:,1:],lat[:,1:],lon[:,:-1],lat[:,:-1])
    dy[1:,:]=haversine(lon[1:,:],lat[1:,:],lon[:-1,:],lat[:-1,:])

    dx[:,0] = dx[:,1]
    dy[0,:] = dy[1,:]
    
    dx = dx * 10**3
    dy = dy * 10**3

    area = dx*dy
    grid_distance = np.mean(np.append(dy[:, :, None], dx[:, :, None], axis=2))

    return dx,dy,area,grid_distance