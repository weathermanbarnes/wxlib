Basic concepts
==============

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
   much more straightforward if these data structures are available. 
 * File input/output. There are already functions available for python that can read pretty much any
   data format, including netCDF, mat, and GRIB files. There's no need to reinvent those, or even to
   implement your own functions in Fortran.
 * Plot functions.

A note on NCL
"""""""""""""

There are alternative high-level languages, which provide the same advantages over Fortran as python,
and that are similarly easlily coupled with Fortran. The most notable alternative is NCL. In principle, 
it would be possible (and desirable) to provide all Fortran functions in dynlib also in NCL, but 
unfortunately NCL does not support the Fortran90 concept of modules. Until this problem is solved, 
dynlib cannot be used with NCL.


.. _sec_examples:

Examples: Read data, apply diagnostic, plot and save results
------------------------------------------------------------

.. literalinclude:: inc/basic_example_inst.py
   :linenos:
   
The first two lines are a standard header for a python file. Line 1 declares the file to be a python 
script, allowing it to be executed by ``./script_name.py``. Line 2 gives the file encoding, which is 
required of you use non-acsii characters in your script. Both lines are optional, but recommended. 

Lines 4 and 5 import some standard functionality from dynlib. :mod:`dynlib.shorthands` contains some
widely used functions such that they can be imported together and ``conf`` contains some default 
settings. 

The concept behind lines 7 and 8 will be discussed more thoroughly in the following sections. For now,
it is enough to note that these lines ensure that the script uses ERA-Interim reanalysis as its data base,
and that it wants to use variables that are derived through dynlib.

The last import on line 10 finally makes all functions of :mod:`dynlib.diag` available.

The time interval on line 13 is a list of datetime objects, that define the time interval for which
data is to be fetched in lines 17 and 18. Analogously the vertical level for the requested data is
set on line 14. In additon to the actual data, the data fetcher function :func:`get_instantaneous` also 
returns meta-information about the variable, mainly the grid that the variable is defined on.

Line 21 finally calculated total deformation from the given wind velocity components. The results are
then saved to a netCDF file in line 25.

For each of the time steps available, the deformation field is also plotted on a map covering the 
North Atlantic (line 29/30).

.. literalinclude:: inc/basic_example_year.py
   :linenos:
   :emphasize-lines: 13,17-18

If you want to calculate deformation for all the data in an input file, the structure of the script 
is very similar. The differences are highlighted in yellow. Most notably, instead of using
:func:`get_instantaneous` the script calls :func:`metopen` which simply returns all the data available 
for the requested variable in the requested file. 

File names typically follow a standard structure, which is stored in the configuration key 
``conf.file_std``. Parameters for completing the file name are filled in using the python string
formatting mechanism with the ``%`` operator. This mechanism might seem a bit tedious at first, but it
will make sure that your scripts run unchanged, when for example the file naming convention changes
or when you apply the script to a different data set.


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


Overview over available contexts
""""""""""""""""""""""""""""""""

Currently, the following contexts exist:

:mod:`dynlib.context.erainterim`
   ERA-Interim Reanalysis, 0.5° resolution, covering 1979-almost present.
:mod:`dynlib.context.ecmwf_fc`
   ECWMF forecast archive for the last 2 weeks or so. Currently contains only the deterministic forecast at 0.5° resolution.
:mod:`dynlib.context.metno_fc`
   AROME2.5 forecast archive for the last 2 weeks or so. 2.5 km resolution covering most of Scandinavia.
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
:mod:`dynlib.context.nora10`
   NORA10 local reanalysis for Norway, 10 km resolution.
:mod:`dynlib.context.bedymo`
   Bedymo model output.
:mod:`dynlib.context.wrf`
   WRF model output.


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

