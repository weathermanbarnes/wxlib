Version 0.30.0
==============

Key additions
-------------

 * Python 3 support

 * Context switches

 * Some new diagnostics and utilities, for example:
   
   - Filling NaN-values (-> Sebastian)
   - Labelling connected regions to define objects
   - Potential temperature from temperature and pressure (-> Lars)


Version 0.20.0 "Newsettings"
============================

Key backwards incompatible changes
----------------------------------

 * New structure for plot settings. New syntax and, the settings can now differ 
   between vertical levels.

   Old version
      conf.contourf.<q>[<key>] = <value>
   New version
      conf.plotf[<plev>,<q>,<key>] = <value>

   plev and q can be None to set settings for all vertical levels/all variables 
   at a certain level or access this default setting.

   Analogously conf.contour transitioned to conf.plot.

 * Changed filename pattern specification

   Old version
      conf.file_std % {'time': 2011, 'plev': 'pv2000', 'q': 'u'}
   New version
      conf.file_std % {'time': 2011, 'plev': 'pv2000', 'qf': 'u'}

   For consistency conf.qi is now called conf.qf. The 'time' keyword does now
   expect strings. Integers or other types are still allowed, and implicitly
   converted to string.

   conf.file_stat and conf.file_mstat do not exist anymore. conf.file_mstat is 
   now superseeded by the more general pattern conf.file_agg, which is however 
   not yet fully incorporated in all parts of dynlib.

   Similarly, for future use the file patterns conf.file_comp for composites, 
   conf.file_eof for EOFs and conf.file_ts for time series have been added.

 * Changed default configuration

   Previously, the default configuration was very specific to ERA-Interim 
   reanalysis. The new default settings do not contain any data set specific 
   information.
   
   The previous default settings can be restored by using the erainterim, 
   derived and plot contexts::

      >>> import dynlib.context.erainterim
      >>> import dynlib.context.derived
      >>> import dynlib.context.plot

   Refer to the documentation for more information about the new concept of 
   contexts.

 * Renamed overlays:
   
   map_overlay_dat is now named map_overlay_contours and analogously, 
   section_overlay_dat now section_overlay_contours.

 * More general metopen and get_instantaneous

   These functions stay largely compatible with previous versions. The only
   backwards incompatible change is that they no longer support cutting the 
   plot domain in the x- and y-dimensions. This feature was unused, and most
   likely broken since some time.

   metopen has been generalised to allow usage without specifying a variable.
   In this case, no data is returned, but only the file object and potentially
   the static grid information.

   get_instantaneous has been generalised to work with arbitrary input data
   sets, as long as some pertinent meta information is available:

    * conf.timestep (time step in the data set)
    * conf.gridsize (grid size of the data set)
    * conf.times or conf.years (valid keys in the time dimension)

   Furthermore, the time information in the static grid object returned by
   get_instantaneous is now correct and can directly be used for saving in
   metsave.

 * Easier-to-use metsave

   Metsave does no longer take the year/time information, but rather relies
   on the time information within the static grid object. metsave does now
   also respect the name of the x- and y-dimensions as given in the grid
   object. 

 * Proper handling of missing values

   Missing values are now represented by NaNs. A new constant
   dynfor.consts.nan is available to be able to set NaNs explicitly in Fortran
   code. NaNs are correctly transferred between python and Fortran. 

   Derivative functions do now return NaN instead of zero at boundaries where a
   centered difference estimation of the derivative is impossible.

   If data is read/written to a compressed netCDF file with 16-bit integers,
   the scale/unscale functions transparently convert the integer missing
   values to NaN and back.
 
 * Removed obsolete streamplot
   
   This functionality is now provided by matplotlib itself.


Key additions to dynlib
-----------------------

 * Contexts

   Contexts are a modular set of settings, either to describe a certain data 
   set, or for a specific task like plotting. A list of available contexts and 
   a list of planned contexts is included in the documentation.

 * Documentation

   The API documentation is complete, and the general introduction became 
   pretty comprehensive. 

 * New Fortran modules

   All Fortran subroutines that calculate terms in the tendency equations for 
   all kinds of variables are now collected in a new tend module.

 * New set of diagnostics

   Lukas Patritz's slope tool calculating slopes of geopotential on isentropic 
   levels and the different terms in its tendency equation are included in 
   dynlib.


In addition there are many, many smaller changes. Contact me, in case you
experience problems that cannot be explained by the above changes.
