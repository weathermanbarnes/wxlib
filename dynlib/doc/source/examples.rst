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

Lines 4 through 11 import some standard functionality from dynlib. Of these line 6 is of special 
importance, as this is line defines that the script uses ERA-Interim reanalysis as its data base.

The time interval on line 13 is a list of datetime objects, that define the time interval for which
data is to be fetched in lines 17. Analogously the vertical level for the requested data is
set on line 14. In additon to the actual data, the data fetcher function :func:`get_instantaneous` also 
returns meta-information about the variable, mainly the grid that the variable is defined on.

Line 20 finally calculated total deformation from the given wind velocity components. The results are
then saved to a netCDF file in line 26.

For each of the time steps available, the deformation field is also plotted on a map covering the 
North Atlantic (line 30-32).

.. literalinclude:: inc/example_diag_yearly.py
   :linenos:
   :emphasize-lines: 11-12,29

The script changes hardly at all if you want to calculate and save deformation for longer period.
Typically, for ERA-Interim, we would use yearly chunks as shown in the script. Differences to the first
example are highlighted in yellow.  


Define and calculate composites
-------------------------------

.. literalinclude:: inc/example_composite.py
   :linenos:
   :emphasize-lines: 32
   
Requesting composite averages works very similarly to requesting instantaneous data. The only additional
argument in :func:`get_aggregate` over :func:`get_instantaneous` is the list of composite conditions to 
be applied. 

These conditions can be constructed easily from the primitives defined in :mod:`dynlib.metio.composite`. 
In this example, lines 18-19 defines two conditions based on the date only to define seasonal averages 
through a composite. In the following, line 22 loads a monthly NAO index time series, and lines 23-24 
define composites based on the NAO index. Finally, line 28 combines the seasonal and the NAO criteria
by creating all combinations, here composites of NAO+ and NAO- for both summer and winter.

The function  :func:`get_instantaneous` does (currently) not return the grid information itself, so we
need to manually request this metadata through :func:`get_static` (line 35) before we can save the result.

.. literalinclude:: inc/example_composite_datadriven.py
   :linenos:
   :emphasize-lines: 21, 24

The final example adapts the previous composite definitions to an example of a data-driven composite.
This example composites all time steps in which the ERA-Interim grid point closest to Bergen exceeds
a 2-meter temperature of 20°C (line 21). The structure of the script is almost entirely identical to 
the above, only in line 24 the combination of the one seasonal criterion with the one temperature
criterion is done manually instead through the :func:`matrix` function.


