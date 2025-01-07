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
from diag_blocking_indices import blocking_index_masato
from grid import grid_from_xarray

def detect_blocks_by_masato(z500,init_df=None,init_obj_arr=None,
                                  min_overlap_perc=10,min_area=1e6*1e6,min_living_tstep=None,dT=6,
                                  speedcoord="mass_center_coords",
                                  return_data=True,quiet=False,outpath=None,outfilename=None):

    start = time.perf_counter()

    print('Detecting blocks ...')
    print('   - Calculate masato blocking index...')
    b_index=blocking_index_masato(z500)
    grid=grid_from_xarray(b_index)
    
    print('   - Detecting potential blocks (thresh_id)...')
    block_objects = thresh_id(b_index,init_obj_arr=init_obj_arr,
                             threshold=0.,thresh_type='lesser',
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
        
    block_obj_xr=xr.DataArray(
        block_objects,
        dims=dims,
        coords=coords,
        name="blocks",
        attrs={
            'units': 'flag IDs',
            'method': 'masato',
            'description': 'IDs of detected blocks from z500 data by the masato method',
            'minimum area [m2]': min_area_attr,
            'minimum overlap [%]': min_overlap_perc,
            'minimum lifespan [hours]': min_lifespan_attr
        })

    b_index_xr=xr.DataArray(
        b_index,
        dims=dims,
        coords=coords,
        name="blocking_index",
        attrs={
            'units': 'geopotential',
            'method': 'masato',
            'description': 'Blocking index as defined by Masato et al. (2012) from z500',
        })
    
    print('   - Producing the object characteristics ...')
    grACs = calc_object_characteristics(block_obj_xr, b_index, grid)

    if min_area is not None:
        print('   - Remove small highs (by area) ...')
        grACs,block_obj_xr=remove_obj_by_min_area(grACs,block_obj_xr,min_area=min_area)
    
    print('   - Tracking by overlap ...')
    block_df=track_by_overlap(block_obj_xr.values,grACs,grid,min_overlap_perc=min_overlap_perc,
                                 speedcoord=speedcoord,init_obj_arr=init_obj_arr,init_df=init_df)

    if min_living_tstep is not None:
        print('   - Remove short tracks ...')
        block_df,block_obj_xr=remove_short_tracks(block_df,block_obj_xr,timesteps=min_living_tstep,include_final_timestep=False)

    block_dataset = xr.Dataset({
                        "blocks": block_obj_xr,
                        "blocking_index": b_index_xr})
    if outfilename is not None:
        block_dataset.to_netcdf(outpath+outfilename+'.nc',
                                  encoding={"blocks": {"zlib": True, "complevel": 9},
                                            "blocking_index": {"zlib": True, "complevel": 9}})
    
    if outfilename is not None:
        with open(outpath+outfilename+'.pkl', 'wb') as handle:
            pickle.dump(block_df, handle)
    
    print('Completed!')
    end = time.perf_counter()
    timer(start, end)
    
    if return_data:
        return block_obj_xr,block_df,grid,b_index