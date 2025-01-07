import sys
import glob
import numpy as np
import xarray as xr
import pandas as pd
from tqdm import tqdm

from utility import keep_obj_by_ids

def get_clim_count(filenames,invar):
    first=True
    for fn in tqdm(filenames):
        data=xr.open_dataset(fn)[invar]
        data=data.where(data==0,1,0)
        data=data.sum(dim='time')
        if first:
            data_full=data.copy()
            first=False
        else:
            data_full=data_full+data

    return data_full

def get_starting_position_clim_count(filenames_nc,filenames_df,invar):
    df=pd.concat([pd.read_pickle(fn) for fn in sorted(filenames_df)])
    
    track_ids=df.track_id.unique()
    init_df=[]
    for tid in track_ids:
        df_tid=df[df.track_id==tid]
        init_df.append(df_tid[df_tid.time==df_tid.time.min()[0]])
    init_df=pd.concat(init_df)
    
    first=True
    for fn_nc in tqdm(filenames_nc):
        data=xr.open_dataset(fn_nc)[invar]
        init_df_dt=init_df[init_df.time>=data.time[0].values]
        init_df_dt=init_df_dt[init_df_dt.time<=data.time[-1].values]
        data=keep_obj_by_ids(data,init_df_dt.obj_id.to_list())
        data=data.where(data==0,1,0)
        data=data.sum(dim='time')
        if first:
            data_full=data.copy()
            first=False
        else:
            data_full=data_full+data

    return data_full