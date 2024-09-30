Extending dynlib
================

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

Which language to choose for your development?
""""""""""""""""""""""""""""""""""""""""""""""

The discussion above will in most cases give a clear indication on which 
language you should choose for your extension of dynlib. If still in doubt, a good strategy
might be to start implementing in python. When you notice some bottleneck that slows down
the python code too much for your liking, add a specialised Fortran function for that specific
problem.

If you have found the appropriate language, follow the recipies in the following sections to

 * :ref:`sec_add_f_fun`
 * :ref:`sec_add_py_fun`
 * :ref:`sec_add_f_mod`
 * :ref:`sec_add_py_mod`
 * :ref:`sec_add_datasource`

While developing, please make sure to follow the :ref:`sec_f_fmt` and the :ref:`sec_py_fmt`. 
The recipies also explain how to include documentation for your addtions. In case that does not
help your problem, a few further remarks about its general workings are in :ref:`sec_devel_docs`.



Version control system
----------------------

The source code of dynlib is managed using the version control system `git <https://git-scm.com/>`_. 
git is a distributed version control system, which means that technically every copy ("clone in 
git lingo) is equal and hence contains the entire history of the project. There is a 
`git book <https://git-scm.com/book/en/v2>`_ which contains a gentle introduction to git. 
Nevertheless, git will be a bit confusing at first, but it is worth the learning effort!

The following is a quick overview over the most important git commands that will enable you to 
start developing for dynlib. It is not a substitute for the git book!

First you will need to create a clone of the repository, which will be come your personal 
working copy of the source code.::

  $ git clone https://git.app.uib.no/Clemens.Spensberger/dynlib.git

When you have a local copy, make sure everything works as expected by compiling dynlib.::
  
  $ cd dynlib
  $ ./compile

Once you have made some changes, first make sure everything works by compiling again. Then you 
can check what changes have accumulated since the last commit by using ::

  $ git status

or in more detail by using ::

  $ git diff
  $ git difftool -t <diff-tool-of-your-choice>

The latter commands will show your changes to the individual files. The ``difftool`` variant allows
to specify an own command for a comparison program. A nice and simple GUI program available at UiB
is ``meld``. 

You can always discard any uncommited changes by ::

  $ git checkout <file(s)>

Once you are happy with your changes, you want to commit a new version of the source code. To do this
first select which of the changed or new files to include in the new version.::

  $ git add <file(s)>

Finally, you can ::
  
  $ git commit 

to save your changes as a new version. Make sure to include a brief but illustrative description of your
changes. As a short cut, if you want to automatically add all changed (but no new!) files to a commit, 
you can use::

  $ git commit -a


Using version control to collaborate
""""""""""""""""""""""""""""""""""""

At some point you will want to make some of your additions available to other users, or take over
some of their changes to dynlib. The recommended way to do so is to fork the `dynlib repository on
the UiB gitlab <https://git.app.uib.no/Clemens.Spensberger/dynlib>`_, commit your changes to your fork of the project and 
then to send a pull request from there.



Adding functions to dynlib
--------------------------

.. _sec_add_f_fun:

Adding a Fortran subroutine
"""""""""""""""""""""""""""

Here's a brief recipe of how to add a fortran function or subroutine to an existing module.

 #. Add the source code of your new function and subroutine to the module cotains section.
    If the new routine belongs thematically to existing routines, put it somewhere close.
 #. Prepend a thorough documentation of the routine.
 #. Fortran functions are unfortunately not discovered automatically by the automatically
    created API documentation. Hence, add a ``autofunction`` line to the appropriate API
    documentation page in ``doc/source/api/``.

Once dyncal is available, instructions of how the new subroutine can be made usable directly
with dyncal will be included here.


.. _sec_add_py_fun:

Adding a python function
""""""""""""""""""""""""

Here's a brief recipe of how to add a python function to an existing module.

 #. Add the source code of your new function and subroutine to the appropriate python module. 
    If the new routine belongs thematically to existing routines, put it somewhere close.
 #. Include a through documentation as the functions docstring. The function and its docstring
    will be automatically included in the API documentation.

Once dyncal is available, instructions of how the new function can be made usable directly
with dyncal will be included here.


.. _sec_add_f_mod:

Adding a Fortran module
"""""""""""""""""""""""

If your new subroutine does not fit into any of the existing Fortran modules, you can create a new
one. The procecures for this are a bit more involved, but still not difficult.

 #. Take example from one of the existing Fortran modules and create an empty module, 
    including a general description of the contents of the module.
 #. Follow the steps in :ref:`sec_add_f_fun` to add your new subroutine.
 #. Follow the steps in :ref:`sec_add_py_mod` to add a new python module containing the contents
    of your new Fortran module. Take example of one of the other python modules named after
    a Fortran module. It is important to take over and adapt the call of :func:`dynlib.docutil.takeover`
    to match your new Fortran module. This command ensures that your documentation in the Fortran
    source code is made available in python and hence for the automatically generated documentation.
    Furthermore, this python module allows to import your new module by::

       import dynlib.my_new_fortran_module

 #. Take example from one of the existing module documentation pages in ``doc/source/api`` and
    create an API documentation page for your new module. 
 #. Add your new documentation page to ``doc/source/apidoc.rst`` to include it in the overall
    documentation.


.. _sec_add_py_mod:

Adding a python module
""""""""""""""""""""""

If your new function does not fit into any of the existing python modules, you can create a new
one. The procecures for this are a bit more involved, but still not difficult.

 #. Take example from one of the existing python modules and create an empty module, 
    including a general description of the contents of the module.
 #. Follow the steps in :ref:`sec_add_py_fun` to add your new function.
 #. Take example from one of the existing module documentation pages in ``doc/source/api`` and
    create an API documentation page for your new module. 
 #. Add your new documentation page to ``doc/source/apidoc.rst`` to include it in the overall
    documentation.


.. _sec_add_datasource:

Adding a data source
--------------------

Data sources are defined as python submodules of ``dynlib.metio``. Most of the functionality is solved
in general in ``dynlib.metio.datasource``, and only some dataset-specific functions are required, for 
example a function mapping requested data to a list of filenames. Have a look at the existing data sources
and adapt from the source that most closely matches the one you want to add.


.. _sec_f_fmt:

Fortran code conventions
------------------------

The following source code formatting conventions strive to make the code base as 
homogeneous as possible and to reduce errors. None of these conventions is a strict
rule -- just make sure that there's a good reason when breaking them.

 * Use Fortan90 (Fortran free format)
 * Use the kind specification ``ni`` and ``nr`` for all variable definitions and 
   numerals in the code.
 * The bodies of modules, programs, subroutines, functions, etc. should be indented by 2 spaces.
 * The bodies of loops/if-clauses, etc. should be indented by 3 spaces.
 * Comments should be indented as the surrounding source code.
 * Comment out blank lines.
 * Clearly separate the variable definitions from the body of a function or subroutine.
 * Use lower-case names.
 * Separate composite names by underscore.
 * Prepend a documentation section directly before every module, subroutine and function. 
   The documentation section must be marked by lines beginning with ``!@`` and must follow
   the `numpy conventions for formatting those docstrings <https://sphinxcontrib-napoleon.readthedocs.org/en/latest/index.html#google-vs-numpy>`_.
   These docstrings will be parsed to automatically create the API documentation. 


.. _sec_py_fmt:

python code conventions
-----------------------

Again, these conventions aim to create a homogeneous source code base, which is easily 
accessible for people new to dynlib.

 * Indent by four spaces (make sure your editor does not sliently convert them to tabs). As
   indention is part of the python syntax it's important not to mix indention patterns!
 * Use lower-case names.
 * Separate composite names by underscore.
 * Comments should be indented as the surrounding source code.
 * Use python docstrings for documentation, and follow the `numpy conventions for
   formatting those docstrings <https://sphinxcontrib-napoleon.readthedocs.org/en/latest/index.html#google-vs-numpy>`_.


.. _sec_devel_docs:

About the documentation mechanism
---------------------------------

The documentation is created using `sphinx <http://sphinx-doc.org/>`_, the API documenation using
the sphinx extension `napoleon <https://pypi.python.org/pypi/sphinxcontrib-napoleon>`_. The
documentation is written in reStructuredText and resides in the folder ``doc/source``. 

The API documentation makes use of python docstrings. For the Fortran routines, these docstrings
are injected during import using :func:`dynlib.docutils.takeover`.


Documentation helper API
""""""""""""""""""""""""

.. toctree::
   :maxdepth: 2

   api/docutil

