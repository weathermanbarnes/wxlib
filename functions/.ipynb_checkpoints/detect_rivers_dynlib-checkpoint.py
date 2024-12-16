#!/usr/bin/env python

''' 
    Generic threshold tracking. Originally derived from MOAAPS.
    Updated and structured by: Michael A. Barnes (ARC CoE for 21st Century Weather, Monash University)

'''

import numpy as np
import os
from scipy import ndimage
from scipy.spatial import ConvexHull
from utility import clean_up_objects, BreakupObjects, ConnectLon_on_timestep, calc_object_characteristics, timer
from utility import minimum_bounding_rectangle, DistanceCoord
from detect_threshold import thresh_id
import time
import xarray as xr

def detect_rivers(ivte,ivtn,grid,
                      dT=6,IVTthreshold=500, MinTimeIVT=9,
                      AR_MinLen = 2000, AR_Lat = 20, AR_width_lenght_ratio = 2,
                      output_flag_netcdf=True,
                      outpath='',outfilename='rivers_test'):
    print('Detecting atmospheric rivers ...')
    
    start = time.perf_counter()

    print('   - Calculating IVT ...')
    IVT = (ivte ** 2 + ivtn ** 2) ** 0.5

    print('   - Detcting potential rivers (thresh_id)...')
    IVT_objects = thresh_id(IVT,
                             threshold=IVTthreshold,MinTime=MinTimeIVT,dT=dT,
                             connectLon=True,
                             thresh_type='greater',
                             break_up=True)

    print('   -  Refining rivers ...')
    AR_obj = ar_check(IVT_objects,
                     AR_Lat,
                     AR_width_lenght_ratio,
                     AR_MinLen,
                     grid['longitude2D'],
                     grid['latitude2D'])

    if output_flag_netcdf:
        print('   - Writing out netcdf flag file ...')
        xr.DataArray(
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
                'minimum lifetime [hours]': MinTimeIVT,
                'mimimum length [km]': AR_MinLen,
                'mnimum latitude of centroid': AR_Lat, 
                'mimimum length to width ratio': AR_width_lenght_ratio,
            }
        ).to_netcdf(outpath+outfilename+'.nc')

    print('   - Calculating river characteristics ...')
    grACs = calc_object_characteristics(AR_obj, # feature object file
                     IVT,         # original file used for feature detection
                     outpath+outfilename,
                     grid['time'],                   # timesteps of the data
                     grid['latitude2D'],             # 2D latidudes
                     grid['longitude2D'],            # 2D Longitudes
                     grid['grid_distance'],
                     grid['area'])

    print('Completed!')
    end = time.perf_counter()
    timer(start, end)
    
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
    #                 plt.plot(PointsObj[simplex, 0], PointsObj[simplex, 1], 'k-')
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
                rgiCenter = np.round(ndimage.measurements.center_of_mass(ObjACT[tt,:,:])).astype(int)
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