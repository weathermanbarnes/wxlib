Reading and writing data
========================

Available data sets
-------------------

Dynlib comes with convenience functions for a couple of widely used datasets. Currently the following 
datasets are defined:

 * ERA5
 * ERA-Interim
 * ERA20-C
 * NORA10

For each of these datasets a set of input/output functions are available. These are 
 
 * ``get_instantaneous``: Request instantaneous data for a time period.
 * ``get_time_average``: Request time-average data for a time period.
 * ``get_aggregate``: Request a sequence of time averages (e.g. pentad or monthly averages) for a time 
   period.
 * ``get_composite``: Request an average of all time steps fulfilling a given criterion within a time 
   period.
 * ``metsave``: Save data to a file, following conventions as far as possible.
 * ``metsave_composite``: Wrapper for metsave that directly takes the data structure as returned by 
   ``get_composite`` and save the data to a file.

These convenience functions make all datasets available under with the same functions, thus making it
straightforward to apply the same analyses to several of these datasets.

For most cases it will be relatively easy to define new datasets and make the same set of functions 
available for them, because all the tedious tasks are solved generically. Thus, for a new dataset, 
the main remaining task will be to implement some dataset-specific functions, like a mapping of the 
requested data to a set of files in which the data will be found.

The implementation of new datasets will be particularly easy if they follow either the ERA5-paradigm 
(shared with ETH) or the ERA-Interim-paradigm (UiB default) of how the data is organised. In these 
cases the implementations for ERA5 and ERA-Interim will only require some few adaptations to the new 
dataset.

.. note:: Instead of using these convencience functions you can always just read the data through
   alternative libraries such as ``xarray`` or ``netCDF4``. When reading data manually, make sure
   to pass the data as ``numpy.ndarray`` objects to dynlib functions, ``xarray.DataArray`` or 
   ``netCDF4.Variable`` objects will most likely not work and might yield unexpected results.


Variable definitions and other meta data describing dataset
-----------------------------------------------------------

For each data source a small set of meta data is available, defining grid resolution, time step,
available variables and such. They are available through the ``conf`` object in the respective
data source module (e.g. ``dynlib.metio.erainterim.conf``).

Data source attributes
""""""""""""""""""""""

Here is a comprehensive overview over the configuration variables that can be defined for each source.
Note that not all attributes are meaningful for all data sources. Note that some configuration variables 
are managed internally using API functions and should not be edited manually. User-editable attributes
are indicated by the check-mark in the edit-column.

=============================== ======= =============================== ======================= =======================
Name                            Edit	Type				Default                 Description
=============================== ======= =============================== ======================= =======================
calendar                                str                             ``'standard'``          CF-compliant calendar name
datapath	        	✓	list of string			``['.', ]``		List of directory to be searched for data input files.
epoch		        	✓	datetime						Reference time for the data set, typically first available time step. Used for example to anchor aggregation periods.
gridsize	        	✓	2-tuple of int			``(, )``		Grid dimension of the data set.
local_timezone	        	✓	string				``'Europe/Oslo'``	Local time zone identifier. Used for the time information in the changelog of data output files.
oformat		        	✓	string				``'nc'``		In which format to save the output. Currently only supprted: ``'nc'`` for netCDF output.
opath		        	✓	string				``'.'``			Where to save data output files.
q		         		dict				``{}``			Mapping from file name segment to variable name.
qf		         		dict				``{}``			Mapping from variable name to file name segment.
q_avg		         		dict				``{}``			Mapping instantaneous to time-average variable names for those variables that change name when taking a time average.
q_bins		         		dict				``{}``			If binned variable: Mapping from variable name to list of bin boundaries.
q_lines		         		dict				``{}``			If line variable: Mapping from variable name to the associated line offset variable name.
q_long		         		dict				``{}``			Mapping from variable name to long variable name.
q_obj		         		dict				``{}``			If feature mask: Mapping from variable name to the name of the time average without the suffix ``_freq``.
q_std		         		dict				``{}``			Mapping from variable name to the standard (i.e. ECMWF/ERA-Interim) variable name.
q_units		         		dict				``{}``			Mapping from variable name to unit description.
staticfile                              str                             ``''``                  File from which the grid information is to be derived. Ideally a netCDF-file containing the orography for the data source.
timestep	        	✓	timedelta or :class:`tagg`	``[]``			Time step of the data set.
=============================== ======= =============================== ======================= =======================


Adapting data source and variable information your personal settings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You might want to adapt some of these settings, for example specifiying additional folders to look
for data, or specifing the default location for saved netCDF files. You can of course adapt 
``conf.datapath`` and ``conf.opath`` in every script where you use them, but in case you find yourself
adding the same lines to many of your script, follow the below recommended recipe to centralise your
settings.

 #. Create separate personal settings files for a data source (for example called ``my_erai_settings.py``),
    where you import the ``erainterim`` data source and modify the ``conf`` object to your liking. Further,
    use the :func:`dynlib.settings.settings_obj.register_variable` function to define new variables. 

 #. Then you can use 

    >>> from my_erai_settings import *

    to import your adaptation of the ``erainterim`` data source instead of the default.


File formats
------------

Follow these format descriptions to be interoperable with dynlib. For all kinds of data,
the standard format is based on the netCDF standard. 

Line data
"""""""""

The dynlib format for line data in netCDF files is compatible with the netCDF3-classic 
standard. 

The main challenge for storing lines in the netCDF format is the non-gridded nature of 
lines. The coordinates of each point along each line are to be stored at sub-grid scale
resolution, such that storing lines as binary masks is not an option. Furthermore, this
format would waste huge parts of the file with zeros.

Dynlib stores lines as lists of coordinates for every time step. The starting indexes for 
each line is stored in a separate variable, whose name usually ends on ``off``, short for 
offset. Hence, every line variable must be of type float and the dimensions:

 * Time
 * Point index 
 * Information type
   
Information type  must be of length 3. The first information is the xidx, second the 
yidx, and third an arbitrary payload. The indexes follow the Fortran convention and start
with 1.0. The payload depends on the line type. For example, jet axes store the wind speed 
at the jet axis in the third component. 

The line offset variable must be of type integer and have the dimensions:

 * Time
 * Line index

The line offset for time ``tidx`` and line index ``lidx`` gives the point index of the first
point belonging to this line. 

Lines following this data format can be direcly overlayed maps, using the 
:func:`dynlib.figures.map_overlay_lines` function, and can be converted into mask fields with
the :func:`dynlib.utils.mask_lines` function.


Time series
"""""""""""

A time series is a netCDF or equivalent file, that contains a variable ``dates`` and another
variable ``values``. Semantically, the i'th value must be a representative value for the time
period from the i'th to the i+1'th date. Hence, the array ``dates`` must be one item longer 
than the array ``values``. 

To allow several time series with different time axes to be stored in one file, the format
can be extended to the more general names ``PREFIX_dates`` and ``PREFIX_values``. In this case,
dates and values are matched by the ``PREFIX``.


Time-aggregated data
""""""""""""""""""""

Time-aggregated should follow the same conventions as the base data. For ERA-Interim, this is
one file per level and variable. Time-aggregated data files may be concatenated in the time
dimension. When concatenating, the length of the time dimension should not be drastically 
longer than that of an input file. For example, for daily values for up to four years of 6-hourly 
data can safely be concatenated into one file.


Composites
""""""""""

Composite files are stored in the netCDF4 format, using netCDF groups to hierarchically group 
variables. Several closely related composites can be stored in the same variables, for example 
positive and negative phases of a variability pattern, or a sequence of time-lag composites.

The group and variable structure within a composite file must follow::

   /$vertical_level/$variable 
   /$vertical_level/${variable}_hist 
   /$vertical_level/${variable}_std 
   /$vertical_level/${variable}_min
   /$vertical_level/${variable}_max

Here, the supplementary information (histogram, standard deviation, minimum pattern, maximum
pattern) is optional.

Except for the histogram, which includes an additional dimension storing bins, all variables
have the dimensions:

 * Composite identifer
 * y-direction
 * x-direction

The composite identifer should be a variable-length string.



All available IO functions
--------------------------

.. toctree::
   :maxdepth: 2   
   
   api/metio_era5
   api/metio_era5_ml
   api/metio_erai
   api/metio_era20c
   api/metio_nora10


Internal IO API
--------------------------

.. toctree::
   :maxdepth: 2   
   
   api/gridlib
   api/metio
   api/composite
   api/settings


