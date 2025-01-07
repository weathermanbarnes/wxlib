#!/usr/bin/env python

''' 
    Generic threshold tracking. Originally derived from MOAAPS.
    Updated and structured by: Michael A. Barnes (ARC CoE for 21st Century Weather, Monash University)

'''

import numpy as np
import time
import xarray as xr
import pandas as pd
import pickle

from utility import minimum_bounding_rectangle, DistanceCoord, timer
from detect_threshold import thresh_id, calc_object_characteristics, remove_obj_by_min_area
from tracking import track_by_overlap, remove_short_tracks
from gridding import make_polar_grid, convert_data_to_polar

def detect_anticyclones_by_zanoms(z_anom,reg_grid,hemis='sh',init_df=None,init_obj_arr=None,
                                  z_anom_threshold=1000,min_overlap_perc=10,min_area=1e6*1e6,min_living_tstep=20,dT=6,
                                  speedcoord="mass_center_coords",
                                  return_data=True,quiet=False,outpath=None,outfilename=None):

    start = time.perf_counter()

    print('Detecting blocks ...')

    print('   - Regridding data to polar stereographic blocks ...')
    m,grid=make_polar_grid(hemis,reg_grid,ratio_new_grid=0.5)
    z_regrid=convert_data_to_polar(hemis,grid,reg_grid,z_anom)
    
    #start = time.perf_counter()
    
    print('   - Detecting potential blocks (thresh_id)...')
    high_objects = thresh_id(z_regrid,init_obj_arr=init_obj_arr,
                             threshold=z_anom_threshold,thresh_type='greater',
                             connectLon=True, break_up=False, MinTime=None, dT=dT)
    
    if grid['latitude'] is None:# or grid['longitude'] is None:
        dims=['time', 'y', 'x']
        coords={'time': grid['time'],
                'latitude': (('y', 'x'), grid['latitude2D']),  
                'longitude': (('y', 'x'), grid['longitude2D']),}
    else:
        dims=['time', 'latitude', 'longitude']
        coords={
            'time': grid['time'],
            'latitude': grid['latitude'],
            'longitude': grid['longitude'],
        }

    if min_living_tstep is not None:
        min_lifespan_attr=min_living_tstep*dT
    else:
        min_lifespan_attr=0

    if min_area is not None:
        min_area_attr=min_area*1
    else:
        min_area_attr=0
        
    high_obj_xr=xr.DataArray(
        high_objects,
        dims=dims,
        coords=coords,
        name="highs",
        attrs={
            'units': 'flag IDs',
            'description': 'IDs of detected anticyclonic z500 anomalies',
            'z anomaly threshold [gpm]': z_anom_threshold,
            'minimum area [m2]': min_area_attr,
            'minimum overlap [%]': min_overlap_perc,
            'minimum lifespan [hours]': min_lifespan_attr
        })
    
    print('   - Producing the object characteristics ...')
    grACs = calc_object_characteristics(high_obj_xr, z_regrid, grid)

    if min_area is not None:
        print('   - Remove small highs (by area) ...')
        grACs,high_obj_xr=remove_obj_by_min_area(grACs,high_obj_xr,min_area=min_area)
    
    print('   - Tracking by overlap ...')
    high_df=track_by_overlap(high_obj_xr.values,grACs,grid,min_overlap_perc=min_overlap_perc,
                                 speedcoord=speedcoord,init_obj_arr=init_obj_arr,init_df=init_df)

    if min_living_tstep is not None:
        print('   - Remove short tracks ...')
        high_df,high_obj_xr=remove_short_tracks(high_df,high_obj_xr,timesteps=min_living_tstep,include_final_timestep=False)
    
    if outfilename is not None:
        high_obj_xr.to_netcdf(outpath+outfilename+'.nc',
                                  encoding={"highs": {"zlib": True, "complevel": 9}})
    
    if outfilename is not None:
        with open(outpath+outfilename+'.pkl', 'wb') as handle:
            pickle.dump(high_df, handle)
    
    print('Completed!')
    end = time.perf_counter()
    timer(start, end)
    
    if return_data:
        return high_obj_xr,high_df,grid,z_regrid