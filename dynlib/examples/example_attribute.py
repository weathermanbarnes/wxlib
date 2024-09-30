import numpy as np
import xarray as xr
from itertools import product
from dynlib import utils # Assuming detect module is imported from a specific module in your codebase

# Open datasets of features
cyc = xr.open_dataset('/Data/gfi/stormrisk/CESM/cycmask/CESM2_LE.BHISTsmbb.LE2-1111.006.ans.1990.sfc.cyc.nc')
mta = xr.open_dataset('/Data/gfi/stormrisk/CESM/mta/CESM2_LE.BHISTsmbb.LE2-1111.006.ans.1990.sfc.mta.nc')
fid = xr.open_dataset('/Data/gfi/stormrisk/CESM/frovo_id/CESM2_LE.BHISTsmbb.LE2-1111.006.ans.1990.850.frovo_id.nc')
cao = xr.open_dataset('/Data/gfi/stormrisk/CESM/cao/CESM2_LE.BHISTsmbb.LE2-1111.006.ans.1990.sfc.cao.nc').sel(time='1990')
lsp = xr.open_dataarray('/Data/gfi/stormrisk/CESM/attribution/test_totP.nc')


#Make MTAs to masks rather than lines. 
mta_mask = utils.mask_lines_with_data( mta.mta_axis.values, mta.mta_aoff.values.astype(int), (cyc.lat.size, cyc.lon.size))

# Create a DataArray with the same dimensions and coordinates
dims = list(cyc.dims)
coords = {dim: cyc[dim].values for dim in dims}

# Rename varaiables to 
da = xr.Dataset(None, coords = cyc.coords)
da['tp'] = lsp
da['A'] = xr.DataArray(data = mta_mask!=0, dims = dims,coords = coords, name = 'mta')
da['F'] = cyc.cycmask>0
da['C'] = fid.frovo_id>0
da['CAO'] = cao.cao

# Rename dimensions
da = da.rename({'lat':'latitude', 'lon' : 'longitude'})

# Create masks conditions (letter should correspond to the variable in the dataset)
masks_conditions = [
    ('F', lambda x: (x.F > 0)),
    ('C', lambda x: (x.C > 0)),
    ('A', lambda x: (x.A > 0)),
    # Add more masks as needed
    # We don't include CAOs as they are filtered out 
]

# Filter out the CAO precipitation prior to attributing 
_da = (da.CAO>3) & (da.F<0.7)
da['tot_precip'] = da.where(~_da).tp

# Do the attribution (needs to be negtative)
precip_attribution = utils.attribute_to_features(da, masks_conditions, var = 'tot_precip', factor=-1)

# Filtered out the CAO precipitation prior to the watershedding - put it back in 
precip_attribution['CAO'] = da.where((da.CAO>3) & (da.C==0)).tp
precip_attribution['CCAO'] = da.where((da.CAO>3) & (da.C>0)).tp