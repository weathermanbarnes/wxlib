##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
######################################## Generic region by threshold identification ##############################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

import numpy as np
import os
from scipy import ndimage
from utility import clean_up_objects, BreakupObjects, ConnectLon_on_timestep

def thresh_id(data,
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
