""""""
"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Cutoff calculation
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import xarray as xr
import geopandas as gpd
from shapely.geometry import Polygon

from wavebreaking.utils.index_utils import calculate_properties, transform_polygons
from wavebreaking.utils.data_utils import (
    get_dimension_attributes,
    check_argument_types,
    correct_dimension_orientation,
)
from wavebreaking.indices.contour_index import decorator_contour_calculation


@check_argument_types(["data"], [xr.DataArray])
@get_dimension_attributes("data")
@decorator_contour_calculation
def calculate_cutoffs(
    data,
    contour_level,
    contours=None,
    min_exp=5,
    intensity=None,
    periodic_add=120,
    *args,
    **kwargs
):
    """
    Identify cutoff structures.
    Dimension names ("time_name", "lon_name", "lat_name"), size ("ntime", "nlon", "nlat")
    and resolution ("dlon", "dlat") can be passed as key=value argument.
    Before the index calculation, the contour lines are calculated if not provided.

    Parameters
    ----------
        data : xarray.DataArray
            data for the contour and cutoff calculation
        contour_levels : array_like
            levels for contour calculation
        contours : geopandas.GeoDataFrame, optional
            contours calculated with wavebreaking.calculate_contours(...,
            original_coordinates=False)
        min_exp : int or float, optional
            Minimal longitudinal expansion of a cutoff event
        intensity : xarray.DataArray, optional
            data for the intensity calculation (hint: use wb_spatial.calculate_momentum_flux)
        periodic_add: int or float, optional
            number of longitudes in degrees to expand the dataset
            to correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0

    Returns
    -------
        cutoffs: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the cutoffs events:
                * "date": date of the cutoffs
                * "level": level of the contour line
                * "com": center of mass in the format (x,y)
                * "mean_var": mean of the variable used for the contour calculation
                * "event_area": area of a cutoff event
                * "intensity": sum of the intensity (momentum flux)
                * "geometry": (Multi)Polygon with the coordinates in the format (x,y)
    """

    # correct dimension orientation if needed
    data = correct_dimension_orientation(data, *args, **kwargs)

    # filter contours from contour iteration
    if contours is None:
        contours = kwargs["contours"]
    contours = contours[
        (contours.exp_lon < contours.exp_lon.max())
        & (contours.exp_lon >= min_exp)
        & contours.closed
    ].reset_index(drop=True)

    # define Polygons
    polys = [Polygon(row.geometry) for index, row in contours.iterrows()]
    gdf = gpd.GeoDataFrame(contours[["date", "level"]], geometry=polys)
    gdf = gdf.reset_index().rename(columns={"index": "id"})

    # calculate properties and transform polygons
    return gpd.GeoDataFrame(
        calculate_properties(gdf, data, intensity, periodic_add, **kwargs),
        geometry=transform_polygons(gdf, data, **kwargs).geometry,
    )
