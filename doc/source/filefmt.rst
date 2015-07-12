File formats
============

Follow these format descriptions to be interoperable with dynlib. For all kinds of data,
the standard format is based on the netCDF standard. 

Line data
---------

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
``dynlib.figures.map_overlay_lines`` function, and can be converted into mask fields with
the ``dynlib.utils.mask_lines`` function.


Time series
-----------

A time series is a netCDF or equivalent file, that contains a variable ``dates`` and another
variable ``values``. Semantically, the i'th value must be a representative value for the time
period from the i'th to the i+1'th date. Hence, the array ``dates`` must be one item longer 
than the array ``values``. 

To allow several time series with different time axes to be stored in one file, the format
can be extended to the more general names ``PREFIX_dates`` and ``PREFIX_values``. In this case,
dates and values are matched by the ``PREFIX``.


Time-aggregated data
--------------------

Time-aggregated should follow the same conventions as the base data. For ERA-Interim, this is
one file per level and variable. Time-aggregated data files may be concatenated in the time
dimension. When concatenating, the length of the time dimension should not be drastically 
longer than that of an input file. For example, for daily values for up to four years of 6-hourly 
data can safely be concatenated into one file.


Composites
----------

*Note*: The format for composites is not final and may still change!

To reduce the amount of files, composite files are stored in the netCDF4 format. This format
allows to hierarchically group variables. Several closely related composites can be stored in 
the same file, for example positive and negative phases of a variability pattern.

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


EOFs
----

*Note*: The format for EOFs is not final and may still change!

To reduce the amount of files, EOF files are stored in the netCDF4 format. This format
allows to hierarchically group variables. Several closely related composites can be stored in 
the same file, for example positive and negative phases of a variability pattern.

The group and variable structure within a composite file must follow::

   dates
   values
   regressed_dates
   regressed_values
   pattern
   /$vertical_level/$regressed_variable 

Time series follow the (extended) time series convention, otherwise all variables have the dimensions:

 * EOF number
 * y-direction
 * x-direction

The EOF numbering should start with 1. REOFs are to be stored in a different files with an analogous 
structure.

