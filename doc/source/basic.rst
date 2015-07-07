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
is very similar. The differences are highlighted in yellow. Most notably, insteand of using
:func:`get_instantaneous` the script calls :func:`metopen` which simply returns all the data available for 
the requested variable in the requested file. 

File names typically follow a standard structure, which is stored in the configuration key 
``conf.file_std``. Parameters for completing the file name are filled in using the python string
formatting mechanism with the ``%`` operator. This mechanism might seem a bit tedious at first, but it
will make sure that your scripts run unchanged, when for example the file naming convention changes
or when you apply the script to a different data set.


Settings and contexts
---------------------

 * Motivation
 * Cascade: Default settings, context settings, own default settings, script overrides
 * Contexts: Different data sets, different needs

