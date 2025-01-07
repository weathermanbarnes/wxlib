##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
######################################## Generic region by threshold identification ##############################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

import numpy as np
import xarray as xr
import os
from scipy import ndimage
import time
from tqdm import tqdm
from utility import timer, switch_lons_to_360, remove_obj_by_ids

def thresh_id(data,init_obj_arr=None,
                threshold=0, thresh_type='greater',
                connectLon=True,
                break_up=True ,MinTime=4, dT=6, **kwargs):

    ''' Detect coherent regions by a specified threshold (https://github.com/AndreasPrein/MOAAP.git)
    Written by: Michael Barnes (ARC Coe for 21st Century Weather, Monash University)
    Based on : MOAAP (Andreas Prein, https://github.com/AndreasPrein/MOAAP.git)
    
    The algorithm detects coherent regions given a specified threshold. There are options for detecting data
    above or below the threshold as well as using the "break-up" method from MOAAP which breaks up long living 
    objects by extracting the biggest object.
    
    Parameters
    ----------
   
    msl : np.ndarray with shape (nz,ny,nx) and dtype float64
        Any suitable field for regions to be detected in.
    threshold : float64
        Required, default ``0``. Threshold to be used in the detection.
    thresh_type : str of either ``greater`` or ``lesser``
        Required, default ``greater``. Defines whether the identified regions should be less than or greater 
        than the speciied threshold.
    connectLon : bool
        *Optional*, defaut ``True``. If true, applies the ``ConnectLon_on_timestep`` connects objects over the 
        datetine on a time-step by time-step basis. Required if using a global grid.
    break_up : bool
        *Optional*, defaut ``True``. Break up long living objects by extracting the biggest object at each time.
    MinTime : int.
        *Optional*, default ``4``. Only required if breakup=True. Minimum lifetime in data timesteps.
    dT : int.
        *Optional*, default ``6``. Only required if breakup=True. Timestep in hours.
   
    Returns
    -------
    id_objects : np.ndarray with shape (nz,ny,nx) and dtype float64
        Array of identified objects, with each object identified by a unique label.
    '''

    if thresh_type=='greater':
        potDATA = (data > threshold)
    if thresh_type=='lesser':
        potDATA = (data < threshold)

    if MinTime is not None:
        rgiObj_Struct=np.zeros((3,3,3)); rgiObj_Struct[:,:,:]=1
    else: 
        rgiObj_Struct=np.zeros((3,3,3)); rgiObj_Struct[1,:,:]=1
        
    out_objects, nr_objectsUD = ndimage.label(potDATA, structure=rgiObj_Struct)
    
    if init_obj_arr is not None:
        if type(init_obj_arr) is type(xr.DataArray()):
            out_objects[out_objects!=0] += np.amax(init_obj_arr.values)
        else:
            out_objects[out_objects!=0] += np.amax(init_obj_arr)
        
    if MinTime is not None:
        print('        '+str(nr_objectsUD)+' object found')
        out_objects, _ = clean_up_objects(out_objects,
                                       dT,
                                min_tsteps=int(MinTime/dT))
        if break_up:
            print('        break up long living objects that have many elements')
            out_objects, object_split = BreakupObjects(out_objects,
                                         int(MinTime/dT),
                                        dT)
    if connectLon:
        print('        connect objects over date line')
        out_objects = ConnectLon_on_timestep(out_objects)

    return out_objects

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
############################################# Threshold utility functions ###################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

def remove_small_short_objects(objects_id,
                               area_objects,
                               min_area,
                               min_time,
                               DT,
                               objects = None):
    """Checks if the object is large enough during enough time steps
        and removes objects that do not meet this condition
        area_object: array of lists with areas of each objects during their lifetime [objects[tsteps]]
        min_area: minimum area of the object (km2)
        min_time: minimum time with the object large enough (hours)
        DT: time step of input data [hours]
        objects: object slices - speeds up processing if provided
    """

    #create final object array
    sel_objects = np.zeros(objects_id.shape,dtype=int)

    new_obj_id = 1
    for obj,_ in enumerate(area_objects):
        AreaTest = np.nanmax(
            np.convolve(
                np.array(area_objects[obj]) >= min_area * 1000**2,
                np.ones(int(min_time/ DT)),
                mode="valid",
            )
        )
        if (AreaTest == int(min_time/ DT)) & (
            len(area_objects[obj]) >= int(min_time/ DT)
        ):
            if objects == None:
                sel_objects[objects_id == (obj + 1)] =     new_obj_id
                new_obj_id += 1
            else:
                sel_objects[objects[obj]][objects_id[objects[obj]] == (obj + 1)] = new_obj_id
                new_obj_id += 1

    return sel_objects

def clean_up_objects(DATA,
                     dT,
                     min_tsteps = 0,
                     obj_splitmerge = None):
    """ Function to remove objects that are too short lived
        and to numerrate the object from 1...N
    """
    
    object_indices = ndimage.find_objects(DATA)
    MaxOb = np.max(DATA)
    MinLif = int(24 / dT)  # min lifetime of object to be split
    AVmax = 1.5

    id_translate = np.zeros((len(object_indices),2))
    objectsTMP = np.copy(DATA)
    objectsTMP[:] = 0
    ii = 1
    for obj in range(len(object_indices)):
        if object_indices[obj] != None:
            if object_indices[obj][0].stop - object_indices[obj][0].start >= min_tsteps / dT:
                Obj_tmp = np.copy(objectsTMP[object_indices[obj]])
                Obj_tmp[DATA[object_indices[obj]] == obj+1] = ii
                objectsTMP[object_indices[obj]] = Obj_tmp
                id_translate[obj,0] = obj+1
                id_translate[obj,1] = ii
                ii = ii + 1
            else:
                id_translate[obj,0] = obj+1
                id_translate[obj,1] = -1
        else:
            id_translate[obj,0] = obj+1
            id_translate[obj,1] = -1

    # adjust the directory strucutre accordingly
    obj_splitmerge_clean = {}

    if obj_splitmerge != None:
        id_translate = id_translate.astype(int)  
        keys = np.copy(list(obj_splitmerge.keys()))
        for jj in range(len(keys)):
            obj_loc = np.where(int(list(keys)[jj]) == id_translate[:,0])[0][0]
            if id_translate[obj_loc,1] == -1:
                del obj_splitmerge[list(keys)[jj]]

        # loop over objects and relable their indices if nescessary
        obj_splitmerge_clean = {}
        keys = np.copy(list(obj_splitmerge.keys()))
        core_translate = np.isin(id_translate[:,0], keys.astype(int))
        id_translate = id_translate[core_translate,:]
        for jj in range(len(keys)):
            obj_loc = np.where(int(list(keys)[jj]) == id_translate[:,0])[0][0]
            mergsplit = np.array(obj_splitmerge[keys[jj]])
            for kk in range(id_translate.shape[0]):
                mergsplit[np.isin(mergsplit, id_translate[kk,0])] = id_translate[kk,1]
            obj_splitmerge_clean[str(int(id_translate[obj_loc,1]))] = mergsplit
        
    return objectsTMP, obj_splitmerge_clean

def ConnectLon_on_timestep(object_indices):
    
    """ This function connects objects over the date line on a time-step by
        time-step basis, which makes it different from the ConnectLon function.
        This function is needed when long-living objects are first split into
        smaller objects using the BreakupObjects function.
    """
    
    for tt in range(object_indices.shape[0]):
        EDGE = np.append(
            object_indices[tt, :, -1][:, None], object_indices[tt, :, 0][:, None], axis=1
        )
        iEDGE = np.sum(EDGE > 0, axis=1) == 2
        OBJ_Left = EDGE[iEDGE, 0]
        OBJ_Right = EDGE[iEDGE, 1]
        OBJ_joint = np.array(
            [
                OBJ_Left[ii].astype(str) + "_" + OBJ_Right[ii].astype(str)
                for ii,_ in enumerate(OBJ_Left)
            ]
        )
        NotSame = OBJ_Left != OBJ_Right
        OBJ_joint = OBJ_joint[NotSame]
        OBJ_unique = np.unique(OBJ_joint)
        # set the eastern object to the number of the western object in all timesteps
        for obj,_ in enumerate(OBJ_unique):
            ObE = int(OBJ_unique[obj].split("_")[1])
            ObW = int(OBJ_unique[obj].split("_")[0])
            object_indices[tt,object_indices[tt,:] == ObE] = ObW
    return object_indices


### Break up long living objects by extracting the biggest object at each time
def BreakupObjects(
    DATA,  # 3D matrix [time,lat,lon] containing the objects
    min_tsteps,  # minimum lifetime in data timesteps
    dT,# time step in hours
    obj_history = False,  # calculates how object start and end
    ):  

    start = time.perf_counter()

    object_indices = ndimage.find_objects(DATA)
    MaxOb = np.max(DATA)
    MinLif = int(min_tsteps / dT)  # min lifetime of object to be split
    AVmax = 1.5

    obj_structure_2D = np.zeros((3, 3, 3))
    obj_structure_2D[1, :, :] = 1
    rgiObjects2D, nr_objects2D = ndimage.label(DATA, structure=obj_structure_2D)

    rgiObjNrs = np.unique(DATA)[1:]
    TT = np.zeros((MaxOb))
    for obj in range(MaxOb):  
        if object_indices[obj] != None:
            TT[obj] = object_indices[obj][0].stop - object_indices[obj][0].start
    TT = TT[rgiObjNrs-1]
    TT = TT.astype('int')
    # Sel_Obj = rgiObjNrs[TT > MinLif]

    # Average 2D objects in 3D objects?
    Av_2Dob = np.zeros((len(rgiObjNrs)))
    Av_2Dob[:] = np.nan
    ii = 1
    
    object_split = {} # this directory holds information about splitting and merging of objects
    for obj in tqdm(range(len(rgiObjNrs))):
        iOb = rgiObjNrs[obj]
        if TT[obj] <= MinLif:
            # ignore short lived objects
            DATA[DATA == iOb] = 0
            continue
        SelOb = rgiObjNrs[obj] - 1
        DATA_ACT = np.copy(DATA[object_indices[SelOb]])
        rgiObjects2D_ACT = np.copy(rgiObjects2D[object_indices[SelOb]])
        rgiObjects2D_ACT[DATA_ACT != iOb] = 0

        Av_2Dob[obj] = np.mean(
            np.array(
                [
                    len(np.unique(rgiObjects2D_ACT[tt, :, :])) - 1
                    for tt in range(DATA_ACT.shape[0])
                ]
            )
        )

        if Av_2Dob[obj] <= AVmax:
            if obj_history == True:
                # this is a signle[ object
                object_split[str(iOb)] = [0] * TT[obj]
                if object_indices[SelOb][0].start == 0:
                    # object starts when tracking starts
                    object_split[str(iOb)][0] = -1
                if object_indices[SelOb][0].stop == DATA.shape[0]-1:
                    # object stops when tracking stops
                    object_split[str(iOb)][-1] = -1
        else:
            rgiObAct = np.unique(rgiObjects2D_ACT[0, :, :])[1:]
            for tt in range(1, rgiObjects2D_ACT[:, :, :].shape[0]):
                rgiObActCP = list(np.copy(rgiObAct))
                for ob1 in rgiObAct:
                    tt1_obj = list(
                        np.unique(
                            rgiObjects2D_ACT[tt, rgiObjects2D_ACT[tt - 1, :] == ob1]
                        )[1:]
                    )
                    if len(tt1_obj) == 0:
                        # this object ends here
                        rgiObActCP.remove(ob1)
                        continue
                    elif len(tt1_obj) == 1:
                        rgiObjects2D_ACT[
                            tt, rgiObjects2D_ACT[tt, :] == tt1_obj[0]
                        ] = ob1
                    else:
                        VOL = [
                            np.sum(rgiObjects2D_ACT[tt, :] == tt1_obj[jj])
                            for jj,_ in enumerate(tt1_obj)
                        ]
                        rgiObjects2D_ACT[
                            tt, rgiObjects2D_ACT[tt, :] == tt1_obj[np.argmax(VOL)]
                        ] = ob1
                        tt1_obj.remove(tt1_obj[np.argmax(VOL)])
                        rgiObActCP = rgiObActCP + list(tt1_obj)

                # make sure that mergers are assigned the largest object
                for ob2 in rgiObActCP:
                    ttm1_obj = list(
                        np.unique(
                            rgiObjects2D_ACT[tt - 1, rgiObjects2D_ACT[tt, :] == ob2]
                        )[1:]
                    )
                    if len(ttm1_obj) > 1:
                        VOL = [
                            np.sum(rgiObjects2D_ACT[tt - 1, :] == ttm1_obj[jj])
                            for jj,_ in enumerate(ttm1_obj)
                        ]
                        rgiObjects2D_ACT[tt, rgiObjects2D_ACT[tt, :] == ob2] = ttm1_obj[
                            np.argmax(VOL)
                        ]

                # are there new object?
                NewObj = np.unique(rgiObjects2D_ACT[tt, :, :])[1:]
                NewObj = list(np.setdiff1d(NewObj, rgiObAct))
                if len(NewObj) != 0:
                    rgiObActCP = rgiObActCP + NewObj
                rgiObActCP = np.unique(rgiObActCP)
                rgiObAct = np.copy(rgiObActCP)

            rgiObjects2D_ACT[rgiObjects2D_ACT != 0] = np.copy(
                rgiObjects2D_ACT[rgiObjects2D_ACT != 0] + MaxOb
            )
            MaxOb = np.max(DATA)

            # save the new objects to the original object array
            TMP = np.copy(DATA[object_indices[SelOb]])
            TMP[rgiObjects2D_ACT != 0] = rgiObjects2D_ACT[rgiObjects2D_ACT != 0]
            DATA[object_indices[SelOb]] = np.copy(TMP)

            if obj_history == True:
                # ----------------------------------
                # remember how objects start and end
                temp_obj = np.unique(TMP[DATA_ACT[:, :, :] == iOb])
                for ob_ms in range(len(temp_obj)):
                    t1_obj = temp_obj[ob_ms]
                    sel_time = np.where(np.sum((TMP == t1_obj) > 0, axis=(1,2)) > 0)[0]
                    obj_charac = [0] * len(sel_time)
                    for kk in range(len(sel_time)):
                        if sel_time[kk] == 0:
                            # object starts when tracking starts
                            obj_charac[kk] = -1
                        elif sel_time[kk]+1 == TMP.shape[0]:
                            # object ends when tracking ends
                            obj_charac[kk] = -1

                        # check if system starts from splitting
                        t0_ob = TMP[sel_time[kk]-1,:,:][TMP[sel_time[kk],:,:] == t1_obj]
                        unique_t0 = list(np.unique(t0_ob))
                        try:
                            unique_t0.remove(0)
                        except:
                            pass
                        try:
                            unique_t0.remove(t1_obj)
                        except:
                            pass
                        if len(unique_t0) == 0:
                            # object has pure start or continues without interactions
                            continue
                        else:
                            # Object merges with other object
                            obj_charac[kk] = unique_t0[0]

                    # check if object ends by merging
                    if obj_charac[-1] != -1:
                        if sel_time[-1]+1 == TMP.shape[0]:
                            obj_charac[-1] = -1
                        else:
                            t2_ob = TMP[sel_time[-1]+1,:,:][TMP[sel_time[-1],:,:] == t1_obj]
                            unique_t2 = list(np.unique(t2_ob))
                            try:
                                unique_t2.remove(0)
                            except:
                                pass
                            try:
                                unique_t2.remove(t1_obj)
                            except:
                                pass
                            if len(unique_t2) != 0:
                                obj_charac[-1] = unique_t2[0]

                    object_split[str(t1_obj)] = obj_charac

    # clean up object matrix
    if obj_history == True:
        DATA_fin, object_split =    clean_up_objects(DATA,
                                    dT,
                                    min_tsteps, 
                                    obj_splitmerge = object_split)
    else:
        DATA_fin, object_split =    clean_up_objects(DATA,
                                    dT,
                                    min_tsteps)

    end = time.perf_counter()
    timer(start, end)

    return DATA_fin, object_split

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
    if grid['longitude'] is not None:
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
                if grid['longitude'] is not None:
                    if object_indices[iobj][2].start==0 and object_indices[iobj][2].stop==var_objects.shape[2]:
                        object_slice = np.copy(swith_var_objects[object_indices_switch[iobj]])
                        lat_idx_slice  = object_indices_switch[iobj][1]
                        lon_idx_slice  = object_indices_switch[iobj][2]
                        lat_slice = switchLat[lat_idx_slice, lon_idx_slice]
                        lon_slice = switchLon[lat_idx_slice, lon_idx_slice]
                    
                obj_mass_center = np.array([ndimage.center_of_mass(object_slice[tt,:,:]==(iobj+1)) for tt in range(object_slice.shape[0])])                

                centroid=np.array([lon_slice[int(round(obj_mass_center[0][0])),int(round(obj_mass_center[0][1]))],
                                   lat_slice[int(round(obj_mass_center[0][0])),int(round(obj_mass_center[0][1]))]])

                if grid['longitude'] is not None:
                    if centroid[0]>=grid['longitude'][-1]:
                        centroid[0]=centroid[0]-360
                    elif centroid[0]<grid['longitude'][0]:
                        centroid[0]=centroid[0]+360
                
                this_object_charac = {
                    "obj_id": iobj + 1,
                    "mass_center_coords": centroid,
                    "tot": obj_tot,
                    "min": obj_min,
                    "min_coords": obj_min_loc,
                    "max": obj_max,
                    "max_coords": obj_max_loc,
                    "mean": obj_mean,
                    "area_m2": obj_size,
                    "time": obj_times,
                    "track_id": None,
                    "track_len": None,
                }

                try:
                    objects_charac[iobj + 1] = this_object_charac
                except:
                    raise ValueError ("Error asigning properties to final dictionary")


        #if filename_out is not None:
        #    with open(outfilename+'.pkl', 'wb') as handle:
        #        pickle.dump(objects_charac, handle)

        return objects_charac

def remove_obj_by_min_area(grACs,obj_arr,min_area=1):
    remove_objs=[]
    for k in grACs.keys():
        if grACs[k]['area_m2']<min_area:
            remove_objs.append(k)

    grACs = {k: v for k, v in grACs.items() if k not in remove_objs}
    #obj_arr[np.isin(obj_arr, remove_objs)] = 0
    obj_arr = remove_obj_by_ids(obj_arr,remove_objs)

    return grACs,obj_arr

