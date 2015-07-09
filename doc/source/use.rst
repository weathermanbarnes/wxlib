Using dynlib
============

Using the dyncal command
------------------------

The dyncal command is still work-in-progress, and will be documented here once in a usable state.


Using with python
-----------------

Dynlib is designed to work well with python. The example scripts in :ref:`sec_examples` give a good
impression of how to use Fortran routines from python. 

A comprehensive list of available functions is in the :ref:`sec_api`.


Using with NCL
--------------

Unfortunately, the NCL mechanism of making Fortran routines usable from within NCL does not work
with the Fortran90 construct of modules. Until this problem is solved, there is unfortunately no
way to use the Fortran routines of dynlib directly from within NCL. 

If you prefer NCL for plotting, but still want to use dynlib for diagnostics or feature detections, 
you will in a first step  need to save them to a netCDF file and can then in a second step load them
into NCL.


Integrate with own Fortran programs
-----------------------------------

The Fortran routines of dynlib can be used in your Fortran programs. The modules of dynlib are 
available to be included via Fortran's ``use`` statement. Once a dynlib module is used, its
routines can be called.

The Fortran ``mod`` files are residing in ``$prefix/include``, the library itself is residing
in ``$prefix/lib/libdynfor.so``. If ``$prefix`` is in the standard search path for the compiler 
and linker, you only need to include the parameter ``-ldynfor`` to the compilation command.

For the centrally installed dynlib at UiB, however, ``prefix=/Data/gfi/users/local``, which is 
not in the standard search path. Hence, you will also need to include the options

``-I$prefix/include``
   Make the compiler find the Fortran mod files.
``-L$prefix/lib``
   Make the linker find the dynfor shared object.
``-Wl,-R$prefix/lib``
   Make the linker embed the non-standard path into the executable, such that it knows there to look 
   for the shared object at runtime.

to the compiler command to make your program compile and run. 

As a summary, here is a complete minimal example, reading the default value of one of dynlib's configuration 
values.

.. literalinclude:: inc/use_fortran_example.f90
   :language: fortran
   :linenos:

Assuming the above file is called ``example.f90``, you can compile it by

.. literalinclude:: inc/use_fortran_compile.sh
   :language: sh
   :linenos:


Handling missing data
---------------------

Internally in dynlib, missing data is represented by ``NaN``, short for Not-a-Number. There are
many advantages of using ``NaN`` over an custom, numerical value like ``-9999.99``:

 #. Arithmetic operations with a numeric missing value will not return the missing value. As a
    consequence one would need to put a where/if-statement around any arithmeric operations to check
    if any of the operands contains a missing value. This procedure is tedious and error-prone.
    In contrast, any arithmetic operation involving a ``NaN`` will yield ``NaN``, without requiring
    any extra code.
 #. It will *just work* with any application. In contrast, any numerical value could be a valid, 
    meaningful value for some specific application, such that the missing value must be kept application 
    dependent.

There are two downsides, however.

 #. At least in Fortran and python there is no integer-``NaN``. Hence, when compressing data to 16-bit
    integer for efficient storage, one must assign a custom missing value. 
 #. Explicitly assigning/creating a ``NaN`` in Fortran, without creating a floating point exception signal 
    (SIGFPE) during run time is comparatively hard. SIGFPEs are generally considered an error, such that
    debugging is made much harder if they are created on purpose.

Dynlib takes care of these downsides. First, the conversion to integer missing values is done 
transparently in the :func:`dynlib.utils.scale` and :func:`dynlib.utils.unscale` functions when
reading/writing compressed data from/to netCDF files. Second, the Fortran library contains a constant 
``dynfor.consts.nan``, where the ``NaN`` value is created during compile time, such that it does not 
create a SIGFPE during run time. Use this value if you want to explicitly assign a missing value in 
Fortran.


