Basic concepts
==============

.. _sec_f_vs_py:

python versus Fortran
---------------------

The main goal behind the combination of python and Fortran is to minimise the total time required
to finish a specific task, including the time for development, the time for debugging, and the
time for actually running the code. python scripts are much quicker to develop and to debug than 
Fortran programs, but calculations tend to be much faster in Fortran. For that reason, anything
that requires serious computation is implemented in Fortran, and anything else in python. The
tool ``f2py`` compiles the Fortran functions such that they are usable from python.

Here, it is important to note that Fortran only provides the functions to be called from python. 
Hence, there is no Fortran program to be executed, but the execution always happens through python.

As an overview, what is written in Fortran?
"""""""""""""""""""""""""""""""""""""""""""

 * All slightly complex calculations, anything going beyond simple arithmetic functions on fields.

And what is written in python?
""""""""""""""""""""""""""""""

 * Simple arithmetic field operations::
   
     y = a*x + b

   Such simple arithmetic expressions are about equally fast in Fortran an python, even if a 
   any or all of the variables are arrays. Hence, there is no need to write a Fortran function
   encapsulating this kind of operation.
 * Anything that requires more complex data structures. A list or a key-value storage are very hard to 
   implement in Fortran, but included in the python language core. Many algorithms can be implemented
   much more straightforwardly if these data structures are available. 
 * File input/output. There are already functions available for python that can read pretty much any
   data format, including netCDF, mat, and GRIB files. There's no need to reinvent those, or even to
   implement your own functions in Fortran.
 * Plot functions.

A note on NCL
"""""""""""""

There are alternative high-level languages, which provide the same advantages over Fortran as python,
and that are similarly easily coupled with Fortran. The most notable alternative is NCL. In principle, 
it would be possible (and desirable) to provide all Fortran functions in dynlib also in NCL, but 
unfortunately NCL does not support the Fortran90 concept of modules. Until this problem is solved, 
dynlib cannot be used with NCL.


Settings and contexts
---------------------

Dynlib comes with predefined settings for a couple of standard data sets, but aims to be easily
configurable to use most other gridded data. The adaptation to different data sets happens 
in **contexts**.

A context is a set of configuration describing pertinent features of the data set. Some of these
features are:

 * Data file location(s)
 * File naming conventions
 * What variables exist on which levels
 * Standard grid type, size and time step

Not all of these features apply to every data set. They can hence be left empty, at the cost that
some of the functionality of dynlib might not work out-of-the-box. However, in those cases, there
is (or should be!) an alternative way of providing the configuration for your specific application.

To keep the settings systems as flexible as possible, there are several layers of where settings
can be set. Each of these layers can overwrite settings defined before.
 
 #. Global defaults for all data sets
 #. Context settings for a data set
 #. Your personal default settings (potentially dependent on active context)
 #. Overrides within a script using dynlib

There are very few global default settings that apply to all data sets. Most of the relevant 
configuration is done in the contexts describing the different data sets. Context can however not
only be used to describe data sets. For examble the context :mod:`dynlib.context.derived` defines
all variables that can be calculated using dynlib, and provides some meta data for them. This
information is for example used, when saving the diagnostics in the above example scripts. 

Furthermore, contexts can provide defaults for plotting. Some global plotting defaults are
defined in :mod:`dynlib.context.plot`. These defaults apply for full variables plotted in color, 
:mod:`dynlib.context.plot.anomaly` and :mod:`dynlib.context.plot.greyscale` provide the respective
alternatives.


.. _sec_context_list:

Overview over available contexts
""""""""""""""""""""""""""""""""

Currently, the following contexts exist:

:mod:`dynlib.context.erainterim`
   ERA-Interim Reanalysis, 0.5° resolution, covering 1979-almost present.
:mod:`dynlib.context.ecmwf_fc`
   ECWMF forecast archive for the last 2 weeks or so. Currently contains only the deterministic forecast at 0.5° resolution.
:mod:`dynlib.context.metno_fc`
   AROME2.5 forecast archive for the last 2 weeks or so. 2.5 km resolution covering most of Scandinavia.
:mod:`dynlib.context.nora10`
   NORA10 local reanalysis for Norway, 10 km resolution.
:mod:`dynlib.context.derived`
   Variable definitions for all diagnostics and detected features available in dynlib.

Plus a few plotting related contexts:

:mod:`dynlib.context.plot`
   Sensible default settings for color plots of full variables (in contrast to anomaly variables).
:mod:`dynlib.context.plot.greyscale`
   Most defaults as in the standard plot context, but color maps and contours are restricted to greyscale.
:mod:`dynlib.context.plot.anomaly`
   Most defaults as in the standard plot context, but color maps and color bar settings appropriate for 
   plotting anomalies.

In addition, a few contexts are reserved for future use:

:mod:`dynlib.context.era40`
   ERA40 reanalysis.
:mod:`dynlib.context.ncep`
   NCEP/NCAR reanalysis.
:mod:`dynlib.context.c20r`
   20th century reanalysis.
:mod:`dynlib.context.bedymo`
   Bedymo model output.
:mod:`dynlib.context.wrf`
   WRF model output.


Context switching
"""""""""""""""""

Sometimes it can become necessary to operate on data from different data sets. For this case, 
dynlib supports loading several contexts in parallel, and switching between them.

For loading contexts describing additional data sets, they need to be imported in the same way as 
the first context was. During the import a new set of configuration is created, and is made 
/active/. It is the /active/ context that affects how dynlib works. So initially, dynlib will 
operate on the newly imported data set. 

You can switch between different contexts using
:func:`dynlib.settings.set_active_context`.

To query if a given set of configuration is /active/, use
:func:`dynlib.settings.is_active_context`. 

Note that there are special contexts, like :mod:`dynlib.context.derived` or anything plot-related, 
that will not create a new set of configuration, but rather modify the currently active context.


Configuration variables
"""""""""""""""""""""""

Here is a comprehensive overview over the configuration variables that can be defined for each context.
Note that some configuration variables are managed internally using API functions and should not be 
edited manually. Those will be adapted automatically when registering variables and vertical levels
via :func:`dynlib.settings.settings_obj.register_variable`. 

While the plot configuration can of course be edited manually (see :ref:`plot configuration`), the 
:class:`dynlib.settings.plot_settings_dict` object itself should not. The available variable/vertical 
level-combinations are again managed through registering variables.

=============================== ======= =============================== ======================= =======================
Name                            Edit	Type				Default                 Description
=============================== ======= =============================== ======================= =======================
datapath	        	✓	list of string			``['.', ]``			List of directory to be searched for data input files.
epoch		        	✓	datetime						Reference time for the data set, typically first available time step. Used for example to anchor aggregation periods.
file_agg	        	✓	string							File naming convention for time-aggregated input files.
file_static	        	✓	string							File name of the static fields file.
file_std	        	✓	string							File naming convention for standard input files.
file_timeless	        	✓	string							File naming convention for input files without time dimension (e.g. composites, eofs).
file_ts			      	✓	string							File naming convention for time series files.
gridsize	        	✓	2-tuple of int			``(, )``			Grid dimension of the data set.
local_timezone	        	✓	string				``'Europe/Oslo'``		Local time zone identifier. Used for the time information in the changelog of data output files.
mlevs		         		list of string			``[]``			Model levels on which data is available.
oformat		        	✓	string				``'nc'``			In which format to save the output. Currently only supprted: ``'nc'`` for netCDF output.
opath		        	✓	string				``'.'``			Where to save data output files.
plevs		         		list of string			``[]``			Pressure levels on which data is available.
plot				(✓)	:class:`plot_settings_dict`				Line/contour plot configuration. See :ref:`plot configuration` for details.
plotf				(✓)	:class:`plot_settings_dict`				Shading/filled-contour plot configuration. See :ref:`plot configuration` for details.
plotformat			✓	string				``'png'``			In which format to save plots. All graphics extensions supprted by matplotlib can be supplied.
plotpath			✓	string				``'.'``			Where to save plots.
ptlevs		         		list of string			``[]``			Potential temperature levels on which data is available.
pvlevs		         		list of string			``[]``			Potential vorticity levels on which data is available.
q		         		dict				``{}``			Mapping from file name segment to variable name.
qf		         		dict				``{}``			Mapping from variable name to file name segment.
q_bins		         		dict				``{}``			If binned variable: Mapping from variable name to list of bin boundaries.
q_long		         		dict				``{}``			Mapping from variable name to long variable name.
q_units		         		dict				``{}``			Mapping from variable name to unit description.
sfclevs		         		list of string			``[]``			Surface levels on which data is available.
times		        	✓	list of int			``[]``			List of time steps available in the data set. See also: ``years``.
timestep	        	✓	timedelta or :class:`tagg`	``[]``			Time step of the data set.
years		        	✓	list of int			``[]``			List of years available in the data set. See also: ``times``.
zlevs		         		list of string			``[]``			Height levels on which data is available.
=============================== ======= =============================== ======================= =======================


Using contexts in your personal settings
""""""""""""""""""""""""""""""""""""""""

You might want to make your personal settings dependent on which context you are in. There are two 
potential mechanisms to achieve this. The choice on which is preferrable will depend mainly on how 
different you want your settings to be in different contexts.

 #. Create separate personal settings files for different contexts

    This approach will be most useful, if there is little or no overlap in your personal settings 
    between the different contexts. Furthermore, you might also want to move the imports of the 
    predefined contexts (lines 5, 7 and 8 in the above examples) into the personal settings file 
    (here called ``mysettings.py``) and then replace line 5 by::

       from mysettings import conf, proj

    and thereby shorten the file header considerably.
 #. Check which context(s) are active from within a settings file

    This approach will be most useful if there is considerable overlap in your settings between 
    different contexts. :mod:`dynlib.settings` defines the function :func:`ìn_context`, which can
    be used to check if the given contexts is/the given contexts are active.

These two mechanisms can be mixed freely.

Furthermore, you can also define your own contexts, to encapsulate for example all relevant settings
for a specific task or paper. To do so, use the function :func:`dynlib.settings.def_context`. Your 
personal contexts can be used in the same way as the pre-defined ones.

