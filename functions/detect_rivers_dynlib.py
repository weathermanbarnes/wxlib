

def atmospheric_river(mflux_u, mflux_v, grid, thresh_ivt=250,thresh_dist=2000, merthresh=None):
    '''
    Atmospheric river detection after Rutz et. al (2014)
    DOI: 10.1175/MWR-D-13-00168.1

    Atmospheric rivers are defined as continous areas of at least 2000 km length with 
    an atmospheric water vapour transport of 250 kg m-1 s-1 (both these threshold can 
    be changed)

    Instead of the standard 250 kg m-1 s-1 threshold a variable threshold can be used (e.g. the 85 percentile)
    Moreover, compared to Rutz et al. (2014) a meridional transport threshold of 50 kg m-1 s-1 is used.

    Parameters
    ----------
    mflux_u : np.ndarray with shape (ny,nx)
        vertical integral of eastward water vapour flux
    mflux_v : np.ndarray with shape (ny,nx)
        vertical integral of northward water vapour flux
    grid :  gridlib.grid
        Grid information as provided by metopen or get_instantaneous.
    thresh_ivt : float
        threshold used in detecting high ivt areas if another variable is used (like tcw) this threshold should be changed. 
        Default : 250 kg m-1 s-1
    thresh_dist : float
        Minimum length (in km) for the area required to be defined as an atmospheric river. Default: 2000
    merthresh : float 
        Meridional transport threshold. Default : False (but common threshold is 50) 
    
    Returns
    -------
    riv_mask : 
        Mask with 1 at the grid points where the atmospheric rivers are detected, and 0 otherwise
    riv_nr : 
        Mask in which each AR is associated with a different number. 
        
    '''
    from scipy.signal import convolve2d
    from scipy.ndimage import label,generate_binary_structure
    #Define output arrays
    riv_mask = np.zeros(mflux_u.shape) 
    riv_nr   = np.zeros(mflux_v.shape)

    #Thresh_ivt should be either one value or an array of similar size as the ivt array

    #1a. Calculate areas with ivt > threshold
    ivt = (mflux_u**2 + mflux_v**2)**0.5
    ivt_mask = (ivt > thresh_ivt)*1.0

    #1b. Nr. the respective areas
    ivt_areas, num_features = label(ivt_mask) #Nr. Areas

    #Correct for areas across the West/Eastern boundaries
    for y in range(ivt_areas.shape[0]):
        if ivt_areas[y, 0] > 0 and ivt_areas[y, -1] > 0:
            ivt_areas[ivt_areas == ivt_areas[y, -1]] = ivt_areas[y, 0]

    #2. Trace the contour by calculating the nr. of neighbours with ivt > threshold
    ivt_mask_conv = convolve2d(ivt_mask,np.ones([3,3]),mode="same",boundary="wrap")

    #3. Loop over different areas to check if maximum distance > distance threshold
    rivnrtemp = 1
    for idx in range(1,num_features+1): 

        templats = grid.y[(ivt_areas == idx) & (ivt_mask_conv < 8.0)]
        templons = grid.x[(ivt_areas == idx) & (ivt_mask_conv < 8.0)]

        #In step 1b some of the areas have been merged, resulting in that
        #Sometimes an area is not existing anymore
        if(len(templats) == 0):
            continue 

        #Boolean to assign if an area is an atmospheric river or not
        riv = False

        #Check if object straddles the equator
        maxlat = np.nanmax(templats)
        minlat = np.nanmin(templats)
        if((maxlat > 0) & (minlat < 0)):
            continue

        #Check if mean poleward transport is above threshold
        if not type(merthresh) ==type(None): 
            mean_mflux_v = np.nanmean(mflux_v[(ivt_areas == idx)])
            if(np.abs(mean_mflux_v) < merthresh):
                continue

        #Calculate maximum difference in lat and lon (in degrees)
        difflat = np.abs(maxlat - minlat) #np.abs(np.nanmin(templats) - np.nanmax(templats))
        difflon = np.abs(np.nanmin(templons) - np.nanmax(templons))
        if(difflon > 180.0):
            difflon = np.abs(difflon -360.0)
        #First check on maximum distance (Check if maximum possible distance < treshold)
        if(((difflat**2 + difflon**2)**0.5)*158 < thresh_dist):
            continue
        #First check if latitude difference is > threshold (then it is an atmospheric river)
        elif(difflat*111.0 > thresh_dist):
            #print("Difflat > thresh")
            riv = True
        #Else: check great circle distance fo all outermost points of area
        else:
            #print("Checking with great circle distance")
            for idx1 in range(len(templats)):
                for idx2 in range(len(templats)-1,-1,-1):  
                    #Calculate great circle distance
                    if(idx1 != idx2):
                        dist = utils.dist_sphere(templats[idx1], templons[idx1], templats[idx2], templons[idx2])
                        #If distance is > threshold distance, stop checking
                        if(dist > thresh_dist):
                            riv= True 
                            break
                if(riv == True):
                    break
        #If the area is an atmospheric river, save to output
        if(riv == True):
            riv_mask[(ivt_areas == idx)] = 1.0
            riv_nr[(ivt_areas == idx)]   = rivnrtemp
            rivnrtemp += 1

    #Return: mask with areas belonging to an atmospheric river (riv_mask)
    #And for every river area a specific nr. 
    return riv_mask, riv_nr