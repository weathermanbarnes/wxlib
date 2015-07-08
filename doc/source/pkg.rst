Packaging and installation
==========================

Software requirements
---------------------

Dynlib uses the following, freely available software.

 * ``git`` (only for version control, not required for installing a packaged version)
 * ``python``
 * ``numpy`` (python package)
 * ``scipy`` (python package)
 * ``matplotlib`` (python package)
 * ``matplotlib-basemap`` (python package)
 * ``Scientific`` (python package for netCDF4 support)
 * ``gfortran`` (Fortran compiler)
 * ``lapack`` (Linear algebra library)
 * ``spherepack`` (Spherical harmonics library)

Spherepack is generally not available as a package and must be compiled and installed manually.
You can download it from the `UCAR website <https://www2.cisl.ucar.edu/resources/legacy/spherepack>`_.


Packaging mechanism: setuptools
-------------------------------

Dynlib follows the standard python packaging and installing mechanism provided by ``setuptools``.
Packaged versions of dynlib can hence be installed with ``pip``. Packages can be created using::

   python setup.by sdist

for a source code packages, without any compiled libraries, and::

   python setup.py bdist

to create a binary package, in which the Fortran library is already compiled. Please note,
however, that binary packages are specific to an operating system.


Install on your own Linux machine
---------------------------------

The general procedure to install dynlib on a Linux machine is:

 #. Install the required software packages using the package manager of your choice. Be aware, though
    that the packages might be named slightly different in different Linux distributions.
 #. Install spherepack manually. 
 #. If you installed spherepack to a location that is not in the default search path for libraries,
    you might need to adapt the ``LIBS=`` definition in the compile script.
 #. Execute ``./compile`` script to ensure the library compiles.
 #. Use ``python setup.py install`` or ``python setup.py install --prefix=/custom/path`` to install
    dynlib.

Please give feedback, in case you encouter trouble, or in case these instructions seem incorrect/outdated.


Install on your own MacOS X machine
-----------------------------------

The overall procedure is identical to Linux, but you have several package managers to choose from.
Dynlib is known to work with the packages installed by `Homebrew <http://brew.sh/>`_, but it should 
also work with `MacPorts <https://www.macports.org/>`_ and others. Be aware, that MacPorts and others
tend to install libraries to a non-standard location, so you might need to adapt paths in dynlib's 
compile script.


Install on your own Windows machine
-----------------------------------

The required software is also available for Windows, so in principle there's nothing in the
way of installing dynlib on a Windows machine. If you succeed in installing dynlib on Windows
please give feedback and your instructions will be included here.


Updating the centrally installed dynlib
---------------------------------------

If you have the appropriate access rights, you can just do::

   python setup.py install --prefix=/Data/gfi/users/local

in your working copy of the dynlib repository to install your version. 

**But even if you might have the access rights: Be careful, and notify other users of 
dynlib of your changes!**

