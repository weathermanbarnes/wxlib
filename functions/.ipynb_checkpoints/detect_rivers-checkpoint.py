#!/usr/bin/env python

''' 
    Generic threshold tracking. Originally derived from MOAAPS.
    Updated and structured by: Michael A. Barnes (ARC CoE for 21st Century Weather, Monash University)

'''

import numpy as np
import os
from scipy import ndimage
from scipy.spatial import ConvexHull
from utility import clean_up_objects, BreakupObjects, ConnectLon_on_timestep, timer #calc_object_characteristics, timer
from utility import minimum_bounding_rectangle, DistanceCoord
from detect_threshold import thresh_id
import time
import xarray as xr
import pandas as pd
import pickle

def detect_rivers(ivte,ivtn,grid,
                      dT=6,IVTthreshold=500,
                      AR_MinLen = 2000, AR_Lat = 20, AR_width_lenght_ratio = 2,
                      return_data=True, quiet=False,
                      outpath=None,outfilename=None,
                      min_overlap_perc=10):
    print('Detecting atmospheric rivers ...')
    
    start = time.perf_counter()

    print('   - Calculating IVT ...')
    IVT = (ivte ** 2 + ivtn ** 2) ** 0.5

    print('   - Detcting potential rivers (thresh_id)...')
    IVT_objects = thresh_id(IVT,
                             threshold=IVTthreshold,MinTime=None,
                             dT=dT,
                             connectLon=True,
                             thresh_type='greater')

    print('   - Refining rivers ...')
    AR_obj = ar_check(IVT_objects,
                     AR_Lat,
                     AR_width_lenght_ratio,
                     AR_MinLen,
                     grid['longitude2D'],
                     grid['latitude2D'])

    print('   - Writing out netcdf flag file ...')
    AR_obj_xr=xr.DataArray(
        AR_obj,
        dims=['time', 'latitude', 'longitude'],
        coords={
            'time': grid['time'],
            'latitude': grid['latitude'],
            'longitude': grid['longitude'],
        },
        name="AR_obj",
        attrs={
            'units': 'flag IDs',
            'description': 'IDs of detected atmopsheric rivers',
            'IVT threshold [kg m-1 s-1]': IVTthreshold,
            'mimimum length [km]': AR_MinLen,
            'mnimum latitude of centroid': AR_Lat, 
            'mimimum length to width ratio': AR_width_lenght_ratio,
        })
    if outfilename is not None:
        AR_obj_xr.to_netcdf(outpath+outfilename+'.nc')

    print('   - Producing the object characteristics ...')
    grACs = calc_object_characteristics(AR_obj, IVT, grid)

    print('   - Tracking by overlap ...')
    AR_df=track_by_overlap(AR_obj,grACs,min_overlap_perc=min_overlap_perc)

    if outfilename is not None:
        with open(outpath+outfilename+'.pkl', 'wb') as handle:
            pickle.dump(AR_df, handle)
    
    print('Completed!')
    end = time.perf_counter()
    timer(start, end)

    if return_data:
        return AR_obj_xr,AR_df

def switch_lons_to_360(longitude,latitude,field):
    # Convert longitude to 0-360
    longitude_360 = (longitude + 360) % 360
    
    # Sort longitudes and rearrange the field accordingly
    sorted_indices = np.argsort(longitude_360)
    longitude_sorted = longitude_360[sorted_indices]
    field_sorted = field[:, :, sorted_indices]  # Rearrange field along the longitude axis

    lon2d,lat2D=np.meshgrid(longitude_sorted,latitude)
    return lon2d,lat2D,field_sorted
    
def ar_check(objects_mask,
             AR_Lat,
             AR_width_lenght_ratio,
             AR_MinLen,
             Lon,
             Lat):

    start = time.perf_counter()
    AR_obj = np.copy(objects_mask); AR_obj[:] = 0.
    Objects=ndimage.find_objects(objects_mask.astype(int))

    aa=1
    for ii in range(len(Objects)):
        if Objects[ii] == None:
            continue
        ObjACT = objects_mask[Objects[ii]] == ii+1
        LonObj = np.array(Lon[Objects[ii][1],Objects[ii][2]])
        LatObj = np.array(Lat[Objects[ii][1],Objects[ii][2]])
        # check if object crosses the date line
        if LonObj.max()-LonObj.min() > 359:
            ObjACT = np.roll(ObjACT, int(ObjACT.shape[2]/2), axis=2)

        OBJ_max_len = np.zeros((ObjACT.shape[0]))
        for tt in range(ObjACT.shape[0]):
            PointsObj = np.append(LonObj[ObjACT[tt,:,:]==1][:,None], LatObj[ObjACT[tt,:,:]==1][:,None], axis=1)
            try:
                Hull = ConvexHull(np.array(PointsObj))
            except:
                ObjACT[tt,:,:] = 0
                continue
            XX = []; YY=[]
            for simplex in Hull.simplices:
                XX = XX + [PointsObj[simplex, 0][0]] 
                YY = YY + [PointsObj[simplex, 1][0]]

            points = [[XX[ii],YY[ii]] for ii in range(len(YY))]
            BOX = minimum_bounding_rectangle(np.array(PointsObj))

            DIST = np.zeros((3))
            for rr in range(3):
                DIST[rr] = DistanceCoord(BOX[rr][0],BOX[rr][1],BOX[rr+1][0],BOX[rr+1][1])
            OBJ_max_len[tt] = np.max(DIST)
            if OBJ_max_len[tt] <= AR_MinLen:
                ObjACT[tt,:,:] = 0
            else:
                rgiCenter = np.round(ndimage.center_of_mass(ObjACT[tt,:,:])).astype(int)
                LatCent = LatObj[rgiCenter[0],rgiCenter[1]]
                if np.abs(LatCent) < AR_Lat:
                    ObjACT[tt,:,:] = 0
            # check width to lenght ratio
            if DIST.max()/DIST.min() < AR_width_lenght_ratio:
                ObjACT[tt,:,:] = 0

        if LonObj.max()-LonObj.min() > 359:
            ObjACT = np.roll(ObjACT, -int(ObjACT.shape[2]/2), axis=2)
        ObjACT = ObjACT.astype(int)
        ObjACT[ObjACT!=0] = aa
        ObjACT = ObjACT + AR_obj[Objects[ii]]
        AR_obj[Objects[ii]] = ObjACT
        aa=aa+1
        
    end = time.perf_counter()
    timer(start, end)

    return AR_obj

# Function to detect overlaps and calculate overlap percentages
def detect_overlap_percentage(current,next_step):
    overlaps = {}  # List to store overlapping pairs and percentages

    # Find unique labels in the current and next time steps
    current_ids = np.unique(current[current > 0])  # Exclude background (label 0)
    next_ids = np.unique(next_step[next_step > 0])

    # For each label in the current time step, check for overlaps
    for cid in current_ids:
        mask_current = current == cid  # Mask for the current label

        for nid in next_ids:
            mask_next = next_step == nid  # Mask for the next label

            # Compute overlap area
            overlap_area = np.sum(mask_current & mask_next)

            if overlap_area > 0:
                # Compute the percentage of overlap
                current_area = np.sum(mask_current)
                next_area = np.sum(mask_next)

                overlap_percentage_current = overlap_area / current_area * 100
                overlap_percentage_next = overlap_area / next_area * 100

                overlaps[nid]={
                    "prev_id": cid,
                    "id": nid,
                    "overlap_area": overlap_area,
                    "prev_overlap_percentage": overlap_percentage_current,
                    "overlap_percentage": overlap_percentage_next
                }

    return overlaps

def calc_object_characteristics(
    var_objects,  # feature object file
    var_data,  # original file used for feature detection
    grid,
    outfilename='',  # output file name and locaiton
    min_tsteps=1,       # minimum lifetime in data timesteps
    split_merge = None  # dict containing information of splitting and merging of objects
    ):
    # ========

    num_objects = int(var_objects.max())
    object_indices = ndimage.find_objects(var_objects)
    switchLon,switchLat,swith_var_objects=switch_lons_to_360(grid['longitude'],grid['latitude'],var_objects)
    object_indices_switch = ndimage.find_objects(swith_var_objects)
    
    if num_objects >= 1:
        objects_charac = {}
        print("            Loop over " + str(num_objects) + " objects")
        
        for iobj in range(num_objects):
            if object_indices[iobj] == None:
                continue
                
            object_slice = np.copy(var_objects[object_indices[iobj]])
            data_slice   = np.copy(var_data[object_indices[iobj]])

            time_idx_slice = object_indices[iobj][0]
            lat_idx_slice  = object_indices[iobj][1]
            lon_idx_slice  = object_indices[iobj][2]

            if len(object_slice) >= min_tsteps:

                data_slice[object_slice!=(iobj + 1)] = np.nan
                grid_cell_area_slice = np.tile(grid['area'][lat_idx_slice, lon_idx_slice], (len(data_slice), 1, 1))
                grid_cell_area_slice[object_slice != (iobj + 1)] = np.nan
                lat_slice = grid['latitude2D'][lat_idx_slice, lon_idx_slice]
                lon_slice = grid['longitude2D'][lat_idx_slice, lon_idx_slice]

                #if min(lon_slice[0])==Lon[0][0] and max(lon_slice[-1])==Lon[0][-1]:
                #    print(Lon[0][-1])


                # calculate statistics
                obj_times = grid['time'][time_idx_slice]
                obj_size  = np.nansum(grid_cell_area_slice, axis=(1, 2))
                obj_min = np.nanmin(data_slice, axis=(1, 2))
                obj_min_loc = np.where(data_slice==obj_min)
                obj_min_loc = np.array([lon_slice[obj_min_loc[1],obj_min_loc[2]][0],
                                       lat_slice[obj_min_loc[1],obj_min_loc[2]][0]])
                obj_max = np.nanmax(data_slice, axis=(1, 2))
                obj_max_loc = np.where(data_slice==obj_max)
                obj_max_loc = np.array([lon_slice[obj_max_loc[1],obj_max_loc[2]][0],
                                       lat_slice[obj_max_loc[1],obj_max_loc[2]][0]])
                obj_mean = np.nanmean(data_slice, axis=(1, 2))
                obj_tot = np.nansum(data_slice, axis=(1, 2))


                # Track lat/lon
                if object_indices[iobj][2].start==0 and object_indices[iobj][2].stop==var_objects.shape[2]:
                    object_slice = np.copy(swith_var_objects[object_indices_switch[iobj]])
                    lat_idx_slice  = object_indices_switch[iobj][1]
                    lon_idx_slice  = object_indices_switch[iobj][2]
                    lat_slice = switchLat[lat_idx_slice, lon_idx_slice]
                    lon_slice = switchLon[lat_idx_slice, lon_idx_slice]
                    
                obj_mass_center = np.array([ndimage.center_of_mass(object_slice[tt,:,:]==(iobj+1)) for tt in range(object_slice.shape[0])])                

                centroid=np.array([lon_slice[int(round(obj_mass_center[0][0])),int(round(obj_mass_center[0][1]))],
                                   lat_slice[int(round(obj_mass_center[0][0])),int(round(obj_mass_center[0][1]))]])
                
                if centroid[0]>=grid['longitude'][-1]:
                    centroid[0]=centroid[0]-360
                elif centroid[0]<grid['longitude'][0]:
                    centroid[0]=centroid[0]+360
                
                this_object_charac = {
                    "obj_id": iobj + 1,
                    "mass_center_coords": centroid,
                    "tot": obj_tot,
                    "min": obj_min,
                    "maincoords": obj_min_loc,
                    "max": obj_max,
                    "max_coords": obj_max_loc,
                    "mean": obj_mean,
                    "size": obj_size,
                    "times": obj_times,
                    "track_id": None,
                }

                try:
                    objects_charac[iobj + 1] = this_object_charac
                except:
                    raise ValueError ("Error asigning properties to final dictionary")


        #if filename_out is not None:
        #    with open(outfilename+'.pkl', 'wb') as handle:
        #        pickle.dump(objects_charac, handle)

        return objects_charac

def track_by_overlap(AR_obj,grACs,min_overlap_perc=10):
    track_ID=0
    for z in range(0,AR_obj.shape[0]):
        IDs=np.unique(AR_obj[z].astype(int))[1::]
        if z==0:
            for ID in IDs:
                track_ID=track_ID+1
                grACs[ID]['track_id']=track_ID
        else:
            overlaps=detect_overlap_percentage(AR_obj[z-1],AR_obj[z])
            for ID in overlaps.keys():
                if overlaps[ID]['overlap_percentage']>=min_overlap_perc:
                    grACs[ID]['track_id']=grACs[overlaps[ID]['prev_id']]['track_id']
                else:
                    track_ID=track_ID+1
                    grACs[ID]['track_id']=track_ID
    
            for ID in list(set(list(IDs)) - set(overlaps.keys())):
                track_ID=track_ID+1
                grACs[ID]['track_id']=track_ID
    
    df=pd.DataFrame(grACs)
    df=df.transpose()

    return df