#!/usr/bin/env python

''' 
    Generic threshold tracking. Originally derived from MOAAPS.
    Updated and structured by: Michael A. Barnes (ARC CoE for 21st Century Weather, Monash University)

'''

import numpy as np
import os
from scipy import ndimage
from scipy.spatial import ConvexHull
from utility import minimum_bounding_rectangle, DistanceCoord, timer, remove_obj_by_ids
from detect_threshold import thresh_id
import time
import xarray as xr
import pandas as pd
import pickle

def detect_overlap_percentage(current,next_step,grid):
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
            overlap_area_in_m2 = np.sum(((mask_current & mask_next)*1)*grid['area'])

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
                    "overlap_area_in_m2": overlap_area_in_m2,
                    "prev_overlap_percentage": overlap_percentage_current,
                    "overlap_percentage": overlap_percentage_next
                }

    return overlaps

def get_track_lengths1(df):
    track_id_counts = df.groupby('track_id').size()
    
    # Update 'track_len' based on conditions
    df['track_len'] = df.apply(
        lambda row: track_id_counts[row['track_id']] if row['track_len'] is None else row['track_len'] + track_id_counts[row['track_id']],
        axis=1
    )

    return df

def get_track_lengths(df):
    # Compute the size of each track_id group
    track_id_counts = df.groupby('track_id').size()
    
    # Define a function to update the 'track_len' column for each group
    def update_track_len(group):
        count = track_id_counts[group.name]
        if group['track_len'].isnull().all():
            group['track_len'] = count
        else:
            max_len = group['track_len'].max(skipna=True)  # Get the max of the group, ignoring NaN
            group['track_len'] = max_len + count
        return group
    
    # Apply the function to each group and retain the original DataFrame structure
    df = df.groupby('track_id', group_keys=False).apply(update_track_len)

    return df

def calc_track_speeds(indf,incoord="mass_center_coords"):

    # Define a function to process each group
    def process_group(df, incoord):
        # Extract longitude and latitude from 'max_coords'
        df["templon"] = df[incoord].apply(lambda x: x[0])
        df["templat"] = df[incoord].apply(lambda x: x[1])
        
        # Calculate distance and time difference
        df["distance_km"] = df.apply(
            lambda row: DistanceCoord(
                row["templon"], row["templat"],
                df.shift().loc[row.name, "templon"],
                df.shift().loc[row.name, "templat"]
            ) if row.name > 0 else 0,
            axis=1
        )
        df["distance_m"] = df["distance_km"] * 1000
        
        df["time_diff_sec"] = df["time"].diff().apply(
            lambda x: x[0].total_seconds() if pd.notnull(x) else None
        )
        
        # Calculate speed in km/h and m/s
        df[f"speed_kmh_{incoord}"] = df["distance_km"] / (df["time_diff_sec"] / 3600)
        df[f"speed_ms_{incoord}"] = df["distance_m"] / df["time_diff_sec"]
        
        return df
    
    # Apply the function to each group using groupby
    processed_df = indf.groupby("track_id").apply(lambda group: process_group(group.copy(), incoord))
    
    # Reset the index if needed
    processed_df.reset_index(drop=True, inplace=True)
    
    # Optional: Select only relevant columns
    processed_df = processed_df[["obj_id",f"speed_kmh_{incoord}", f"speed_ms_{incoord}"]]
    
    return pd.merge(indf, processed_df, on="obj_id")


def track_by_overlap(AR_obj,grACs,grid,
                         init_df=None,init_obj_arr=None,
                         min_overlap_perc=10,speedcoord="mass_center_coords"):
    def drop_duplicates_by_area(indict):
        # Step 1: Group by prev_id
        grouped = {}
        for key, value in indict.items():
            prev_id = value['prev_id']
            if prev_id not in grouped:
                grouped[prev_id] = []
            grouped[prev_id].append((key, value))
        
        # Step 2: Find the entry with the largest overlap_area_in_m2 for each prev_id
        filtered = {}
        for prev_id, entries in grouped.items():
            # Find the entry with the maximum overlap_area_in_m2
            max_entry = max(entries, key=lambda x: x[1]['overlap_area_in_m2'])
            filtered[max_entry[0]] = max_entry[1]
    
        return filtered

    if init_df is None:
        track_ID=0
    else:
        track_ID=init_df.track_id.max()
    for z in range(0,AR_obj.shape[0]):
        IDs=np.unique(AR_obj[z].astype(int))[1::]
        if z==0:
            if init_df is None:
                for ID in IDs:
                    track_ID=track_ID+1
                    grACs[ID]['track_id']=track_ID
            else:
                if type(init_obj_arr) is type(xr.DataArray()):
                    init_obj_arr=init_obj_arr.values
                overlaps=detect_overlap_percentage(init_obj_arr[-1],AR_obj[z],grid)
                overlaps=drop_duplicates_by_area(overlaps)
                for ID in overlaps.keys():
                    if overlaps[ID]['overlap_percentage']>=min_overlap_perc and overlaps[ID]['prev_overlap_percentage']>=min_overlap_perc:
                        #grACs[ID]['track_id']=grACs[overlaps[ID]['prev_id']]['track_id']
                        grACs[ID]['track_id']=init_df[init_df.obj_id==overlaps[ID]['prev_id']].track_id.to_list()[0]
                    else:
                        track_ID=track_ID+1
                        grACs[ID]['track_id']=track_ID
        
                for ID in list(set(list(IDs)) - set(overlaps.keys())):
                    track_ID=track_ID+1
                    grACs[ID]['track_id']=track_ID
        else:
            overlaps=detect_overlap_percentage(AR_obj[z-1],AR_obj[z],grid)
            overlaps=drop_duplicates_by_area(overlaps)
            for ID in overlaps.keys():
                if overlaps[ID]['overlap_percentage']>=min_overlap_perc and overlaps[ID]['prev_overlap_percentage']>=min_overlap_perc:
                    grACs[ID]['track_id']=grACs[overlaps[ID]['prev_id']]['track_id']
                else:
                    track_ID=track_ID+1
                    grACs[ID]['track_id']=track_ID
    
            for ID in list(set(list(IDs)) - set(overlaps.keys())):
                track_ID=track_ID+1
                grACs[ID]['track_id']=track_ID
    
    df=pd.DataFrame(grACs)
    df=df.transpose()

    if init_df is not None:
        df=pd.concat((init_df[init_df.time==init_df.time.max().to_list()[0]][df.keys().to_list()],df))
    
    df=calc_track_speeds(df,incoord=speedcoord)
    df=get_track_lengths(df)

    df=df[df.time>=grid['time'][0]]

    return df

def remove_short_tracks(df,obj_arr,timesteps=4,include_final_timestep=False):
    if not include_final_timestep:
        unended_tracks=df[df.time==df.time.max().values[0]].track_id.to_list()
        short_tracks=df[~df['track_id'].isin(unended_tracks)]
    else:
        short_tracks=df.copy()
    short_tracks=short_tracks[short_tracks.track_len<=timesteps]
    
    df=df[~df['track_id'].isin(short_tracks.track_id.unique())]

    high_obj_xr=remove_obj_by_ids(obj_arr,list(short_tracks.obj_id.to_list()))

    return df,high_obj_xr

def update_init_track_lens(init_df, df):
    # Group by 'track_id' in df and get the max 'track_len' for each group
    max_track_len_df = df.groupby('track_id')['track_len'].max().reset_index()
    
    # Merge the max 'track_len' from df into init_df on 'track_id'
    init_df = init_df.merge(max_track_len_df, on='track_id', how='left', suffixes=('', '_df'))
    
    # Update 'track_len' in init_df with the max 'track_len' from df
    init_df['track_len'] = init_df['track_len_df'].combine_first(init_df['track_len'])
    
    # Drop the temporary 'track_len_df' column
    init_df.drop(columns=['track_len_df'], inplace=True)

    return init_df