#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import absolute_import, unicode_literals
import numpy as np 
import math
from scipy.sparse import csr_matrix
from numpy import loadtxt
from datetime import datetime as dt, timedelta as td

##################################################
# Data processing functions
###################################################

def read_file(st_file, nrskip=1):
    str_id   = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[0],dtype=int)
    str_nr   = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[1],dtype=int)
    str_date = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[2],dtype=int)
    str_lat  = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[4])
    str_lon  = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[3])
    
    str_dt = dt_array(str_date)
    
    return str_id, str_nr, str_dt, str_lat, str_lon

def read_file_clim(st_file, nrskip=0):
    str_id   = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[0],dtype=int)
    str_nr   = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[1],dtype=int)
    str_year = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[2],dtype=int)
    str_day = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[3],dtype=int)
    str_lat  = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[5])
    str_lon  = loadtxt(st_file, comments="#", unpack=False,skiprows=nrskip,usecols=[4])
  
    str_dt = []
    for idx in range(len(str_year)):  

        str_dt.append(dt(str_year[idx]-1,12,1,0) + td(days=(str_day[idx]-1)/4))
	
    str_dt = np.array(str_dt)
    
    return str_id, str_nr, str_dt, str_lat, str_lon
    
#Function to convert timestamp (YYYYMMDDHH) to datetime    
def dt_array(str_date):
    #Determine datetime array for the tracks
    str_dt = []
    str_hour = np.zeros(len(str_date))
    str_day = np.zeros(len(str_date))
    str_month = np.zeros(len(str_date))
    str_year = np.zeros(len(str_date))
    
    for idx in range(len(str_date)):
        year = int(str(str_date[idx])[:4])
        month = int(str(str_date[idx])[4:6])
        day   = int(str(str_date[idx])[6:8])
        hour   = int(str(str_date[idx])[8:10])
        str_hour[idx] = hour
        str_day[idx] = day
        str_month[idx] = month
        str_year[idx] = year
        str_dt.append(dt(year,month,day,hour))
        
    return str_dt

##################################################
# General FUNCTIONS
###################################################

def unnest(l,level=1):
    '''
    Function to unnest list
    
    Input: list l
    Output: Unnested list 
    '''
    
    
    l_unnest = [item for sublist in l for item in sublist]
    
    if(level == 1):
        return l_unnest
    else:
        return unnest(l_unnest,level=level-1)

#Functions to get indices in a large vector
#Code from: https://stackoverflow.com/questions/33281957/faster-alternative-to-numpy-where    
def compute_M(data):
    cols = np.arange(data.size)
    return csr_matrix((cols, (data.ravel(), cols)),
                      shape=(data.max() + 1, data.size))

def get_indices_sparse(data):
    M = compute_M(data)
    return [np.unravel_index(row.data, data.shape) for row in M]

##################################################
# FUNCTIONS TO COMPARE DISTANCES
###################################################
def calc_Rossby_radius(lat=45,N=1.3e-2,H=10):
	return N*H/(2*7.29*10**-5*np.sin(lat*np.pi/180))

#From https://gist.github.com/nickjevershed/6480846
def great_circle(lat1, long1, lat2, long2,dist="kilometers"):
    '''
    Calculates great circle distance between lon1,lat1 and lat2,lon2
    
    '''

    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
        
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
    
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
        
    # Compute spherical distance from spherical coordinates.
        
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
    
    cos = (np.sin(phi1)*np.sin(phi2)*np.cos(theta1 - theta2) + 
           np.cos(phi1)*np.cos(phi2))
    cos[np.where(cos>1)] = 1
    arc = np.arccos( cos )

    # Calculate distance in kilometers or meters
    if (dist == "kilometers"): 
       dist = arc*6378.16
    elif(dist == "meters"):
       dist = arc*6378.16*1000.0

    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return dist

#Compute the spatial and distance between two tracks
def compare_trks(x_1,y_1,t_1,x_2,y_2,t_2,timthresh=None):
    len1 = len(x_1)
    len2 = len(x_2)

    timdiff = np.empty([len1,len2])*np.nan

    lt1 = np.outer(y_1,np.ones(len2))
    ln1 = np.outer(x_1,np.ones(len2))
    lt2 = np.outer(np.ones(len1),y_2)
    ln2 = np.outer(np.ones(len1),x_2)

    avelat = (lt1 + lt2)*0.5
    corrfac = np.abs(calc_Rossby_radius(lat=avelat)) #/calc_Rossby_radius(lat=45))
    dist = great_circle(lt1,ln1,lt2,ln2)/corrfac


    for idx1 in range(len1):
        for idx2 in range(len2): #range(len2): #
            dttemp = ((t_1[idx1] - t_2[idx2]).total_seconds())/3600
            timdiff[idx1,idx2] = dttemp
            
    #Calculate time space diff
    if( timthresh == None):
        timspace_diff = None
    else:
        timspace_diff = (dist**2 + timdiff**2/timthresh**2)**(0.5)
    #t1 = np.outer(t_1,np.ones(len2))
    #t2 = np.outer(np.ones(len1),t_2)
    #timdiff = ((t_1[idx1] - t_2[idx2]).total_seconds())/3600
    if( timthresh == None):
        return dist, timdiff
    else:
        return dist, timdiff, timspace_diff
    
def dist_along_track(x,y):
    lentrk = len(x)
    
    if(lentrk <= 1):
        raise ValueError("Vector must have at least two points")
    else:
        avelat = (y[1:] + y[0:(lentrk-1)])*0.5
        corrfac = np.abs(calc_Rossby_radius(lat=avelat)) 
        dist = great_circle(y[1:],x[1:],y[0:(lentrk-1)],x[0:(lentrk-1)])/corrfac
        
    totaldist = np.nansum(dist)
    
    return dist, totaldist


def connect_cyclones(lons1,lats1,times1,lons2,lats2,times2,
                     Options):
        '''
        Function to check if two storm tracks are 'clustered', according to certain distance criteria
        
        Input:
        lons1,  longitude (in degrees) of track 1
        lats1, latitude (in degrees) of track 1
        times1, times of track 1
        lons2,  longitude (in degrees) of track 2
        lats2, latitude (in degrees) of track 2
        times2, times of track 2
        Options, 
        
        In Options:
	distthresh = 1.0 #1. Distance criterium (in Rossby Radii)
        timthresh = 36.0 #2. Time criterium (in hours)
        lngthresh = 1.5 #3. Length overlap criterium (in Rossby Radii)
        timlngthresh = 48.0 #4. Time overlap criterium (in hours)
        '''
    
        conn = 0
        angle = 0
        dt = 0
        dr = 0
        str_contemp1 = np.zeros(len(lons1))
        str_contemp2 = np.zeros(len(lons2))
    
        dists, timdiffs, timspacediff  = compare_trks(lons2,lats2,times2,lons1,lats1,times1,Options["timthresh"]) 

        ###############################################################################
        # Step 1: Selection of points over which a particular storm tracks are connected
        # First it is checked at each point if the local distance and temporal criteria are satisfied, 
        # gives a matrix of m x n, with m length of track 1, and n length track 2
        ###############################################################################
        point_check = (np.abs(timdiffs) <= Options["timthresh"]) & (dists <= Options["distthresh"])
        #Secondly, check which points along track 1 are connected to any (arbitrary) point of track 2
        pntselect = np.nanmax(point_check,axis=0) #*Rossby_45
        test1 = np.nansum(pntselect) 
        
        #Do the same for the other track
        pntselect2 = np.nanmax(point_check,axis=1) #*Rossby_45
        test2 = np.nansum(pntselect2)

        ###############################################################################
        # Step 2: Calculate distance and time over which storm tracks are connected
        # Only makes sense if at least two points are connected along each track
        ###############################################################################
        if((test1 >=2) & (test2 >= 2)):

            #Just select the points which are connected and calculate distances between the points for both tracks
            owndists, owntims, = compare_trks(lons1[pntselect],lats1[pntselect],times1[pntselect],lons1[pntselect],lats1[pntselect],times1[pntselect])
            owndists2, owntims2, = compare_trks(lons2[pntselect2],lats2[pntselect2],times2[pntselect2],lons2[pntselect2],lats2[pntselect2],times2[pntselect2])
            
            if(Options["distmeth"] == "AlongTracksDirect"):
                #Just select the points which are connected and calculate distances between the points for both tracks
                maxdist = (np.nanmax(owndists) + np.nanmax(owndists2))/2.0
            elif(Options["distmeth"] == "AlongTracks"):
                alongdists1, totaldist1 = dist_along_track(lons1[pntselect],lats1[pntselect])
                alongdists2, totaldist2 = dist_along_track(lons2[pntselect2],lats2[pntselect2])
                maxdist = (totaldist1 + totaldist2)/2.0
                
            #maxdist = (np.nanmax(owndists) + np.nanmax(owndists2))/2.0
            maxtime = (np.nanmax(np.abs(owntims)) + np.nanmax(np.abs(owntims2)))/2.0
            
            ratio = (maxtime/(Options["timlngthresh"]))/(maxdist/Options["lngthresh"])
            angle = np.arctan(ratio)
        else:
            maxdist = 0
            maxtime = 0
        
        if((maxtime >= (Options["timlngthresh"])) or (maxdist >= Options["lngthresh"])):
            if(maxtime >= (Options["timlngthresh"])):
                conn += 2

            if(maxdist >= Options["lngthresh"]):
                conn += 1

            str_contemp1[pntselect] = 1.0
            str_contemp2[pntselect2] = 1.0

            angle = angle*180/np.pi

            if(angle == 0):
                print("Zero angle")
                print((maxdist/Options["lngthresh"]))
                print((maxtime/(Options["timlngthresh"])))

            dr = (maxdist/Options["lngthresh"])
            dt = (maxtime/(Options["timlngthresh"]))
                
        return conn, angle, dt, dr, str_contemp1, str_contemp2
    
#Recursive function to find uniquely connected cluster of storms + Type of cluster
def find_cluster_type_dokm(cluster,connTracks,contype="All"):
    #print("CLustering analysis for the following storms:")
    cluster_old = cluster

    #Loop over storms to find connected storms
    for stridx in cluster: 
        conntemp = connTracks.getrow(stridx).data #-1
        #nonzero  = connTracks.getrow(stridx).nonzero()[1]
        nonzero  = connTracks[stridx,::].nonzero()[1]
        
        if(len(conntemp) > 0):    
            if(contype == "All"):
                strmstemp = np.where(conntemp > 0)[0]  #+ stridx
                #typetemp = conntemp[conntemp > 0] 
                cluster = np.append(cluster,nonzero)
            elif(contype == "Bjerknes" and np.nansum((conntemp == 1.0) | (conntemp == 3.0)) > 0):
                strmstemp = nonzero[np.where((conntemp == 1.0) | (conntemp == 3.0))[0]]  #+ stridx
                cluster = np.append(cluster,np.array(strmstemp,dtype=int))
            elif(contype == "Time" and np.nansum(conntemp >= 2.0) > 0):
                strmstemp = nonzero[np.where(conntemp >= 2.0)[0]]  
                cluster = np.append(cluster,np.array(strmstemp,dtype=int))    
            if(contype == "Stagnant" and np.nansum((conntemp == 2.0)) > 0):
                strmstemp = nonzero[np.where(conntemp == 2.0)[0]]
                cluster = np.append(cluster,np.array(strmstemp,dtype=int))     
            
    #Remove duplicate storms
    cluster = np.unique(cluster)

    #Check if all storms are counted??
    if(len(cluster) == len(cluster_old)):
        return cluster_old
    else:
        #connTypes.extend(list(typetemp))
        return find_cluster_type_dokm(cluster,connTracks,contype=contype)
