.. _sec_examples:

Example scripts
===============

Often it is quickest to learn from examples. Here are some scripts illustrating how to solve some common
tasks using dynlib. They also provide a good basis for your own scripts.


Read data, apply diagnostic, plot and save results
--------------------------------------------------

.. literalinclude:: inc/example_diag_plot_instantaneous.py
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

.. literalinclude:: inc/example_diag_yearly.py
   :linenos:
   :emphasize-lines: 12,16-17

If you want to calculate deformation for all the data in an input file, the structure of the script 
is very similar. The differences are highlighted in yellow. Most notably, instead of using
:func:`get_instantaneous` the script calls :func:`metopen` which simply returns all the data available 
for the requested variable in the requested file. 

File names typically follow a standard structure, which is stored in the configuration key 
``conf.file_std``. Parameters for completing the file name are filled in using the python string
formatting mechanism with the ``%`` operator. This mechanism might seem a bit tedious at first, but it
will make sure that your scripts run unchanged, when for example the file naming convention changes
or when you apply the script to a different data set.


Define and calculate composites
-------------------------------

.. literalinclude:: inc/example_composite.py
   :linenos:
   
For creating the composites, you do not need to load any data explicitly. It is enough to define which
variables on which levels you want to composite (line 11), and which compositing criteria you want to 
use (line 15-21). Dynlib comes with a number of predefined compositing criteria for common tasks. These 
include, amongst others, criteria based on time only to create monthly or seasonal averages and several
climate variability indices. A complete list of the predefined criteria is available in 
:mod:`dynlib.composite.tests`.

The call to ``build`` calculates the actual composites. The composites are based on the time 
range defined in ``conf.years``. Finally, ``save`` saves the results in a netCDF file. 

In addition to the predefined compositing criteria in :mod:`dynlib.composite.tests`, dynlib also provides
templates for other common types of criteria in :mod:`dynlib.composite.decider`. For example, a common
pattern of criteria uses a threshold for a position in a meterological field to define whether a time step
should be part of a composite. The second compositing example script illustrates how these template 
criteria might be used.

.. literalinclude:: inc/example_composite_custom.py
   :linenos:
   :emphasize-lines: 17

The overall structure of the script is identical to the example above. The only difference is the 
definition of a new compositing criterion in line 17. Here, the criterion is the exceedence of a
2m-temperature threshold in the ERA-Interim grid point around Bergen. The first argument to
:func:`dat_lowerbound` is user-defined arbitrary name that will be used, for example, to determine
the file name of the resulting results file. The fourth argument contains the indices for the 
grid point corresponding to Bergen, and the final argument contains the threshold to be applied.


