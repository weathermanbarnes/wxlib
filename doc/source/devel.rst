Extending dynlib
================

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

  $ SHARED=/Data/gfi/users/local
  $ git clone $SHARED/src/dynlib.git

If you are not within the UiB network, you might want to use::

  $ git clone <username>@login.uib.no:$SHARED/src/dynlib.git

to clone the repository over ssh. 

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
------------------------------------

At some point you will want to make some of your additions available to other users, or take over
some of their changes to dynlib. This section briefly outlines a recommended procedure for doing so.

First you will need to setup another copy of your repository. Make sure this repository is readable
by others, as this is where you can make your changes available. From within your working copy, do::

  $ git init --bare /path/to/your/personal/public/repository.git
  $ git remote add <name> /path/to/your/personal/public/repository.git

With these commands, you first setup a "bare" repository at the given location. A bare repository contains
all the versioning information, but it cannot be used as a working copy, because no files are checked out.
Bare repositories are only intended for sharing. Second, you setup a "remote" repository in your working 
copy, which is essentially a name for the bare repository that you just created.

Then, if you want to share the current version of your working copy, use ::

  $ git pull <name> <branch>

If you don't now what branch you're working on, it will most likely be ``master``.

The procedures to take someone else's changes is relatively similar. First you will again want to
setup a remote reposity with a name::
  
  $ git remote add <someothername> /path/to/someone/elses/public/repository.git

Then you can get their most recent published version by::

  $ git pull <someothername>

git will then try to combine your current version with the one you pulled. In some easy cases (for 
example of you edited different files / sections of files than your collaborator) git can do the merge
automatically, but more often there will be *conflicts*. A good way to *resolve* your conflicts is to 
use :: 

  $ git mergetool -t <diff3-tool-of-your-choice>

Again, ``meld`` ia good choice for the diff3-tool. For every file in which there is a conflict, this
command will open a 3-way merge window. The resulting, merged file is in the middle, the to opposing
versions left and right. Go through all the difference between the files to come up with a combined 
version. Once you're finished, make sure to commit the merged version::

  $ git commit 

and to potentially share it with your collaborator.


Which language to choose?
-------------------------

The discussion in :ref:`sec_f_vs_py` will in most cases give a clear indication on which 
language you should choose for your extension of dynlib. If still in doubt, a good strategy
might be to start implementing in python. When you notice some bottleneck that slows down
the python code too much for your liking, add a specialised Fortran function for that specific
problem.

If you have found the appropriate language, follow the recipies in the following sections to

 * :ref:`sec_add_f_fun`
 * :ref:`sec_add_py_fun`
 * :ref:`sec_add_f_mod`
 * :ref:`sec_add_py_mod`
 * :ref:`sec_add_context`

While developing, please make sure to follow the :ref:`sec_f_fmt` and the :ref:`sec_py_fmt`. 
The recipies also explain how to include documentation for your addtions. In case that does not
help your problem, a few further remarks about its general working are in :ref:`sec_devel_docs`.



.. _sec_add_f_fun:

Adding a Fortran subroutine
---------------------------

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
------------------------

Here's a brief recipe of how to add a python function to an existing module.

 #. Add the source code of your new function and subroutine to the appropriate python module. 
    If the new routine belongs thematically to existing routines, put it somewhere close.
 #. Include a through documentation as the functions docstring. The function and its docstring
    will be automatically included in the API documentation.

Once dyncal is available, instructions of how the new function can be made usable directly
with dyncal will be included here.


.. _sec_add_f_mod:

Adding a Fortran module
-----------------------

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
----------------------

If your new function does not fit into any of the existing python modules, you can create a new
one. The procecures for this are a bit more involved, but still not difficult.

 #. Take example from one of the existing python modules and create an empty module, 
    including a general description of the contents of the module.
 #. Follow the steps in :ref:`sec_add_py_fun` to add your new function.
 #. Take example from one of the existing module documentation pages in ``doc/source/api`` and
    create an API documentation page for your new module. 
 #. Add your new documentation page to ``doc/source/apidoc.rst`` to include it in the overall
    documentation.


.. _sec_add_context:

Adding a context
----------------

Any settings file can define a new context by calling :func:`dynlib.settings.def_context`. If you
created a new context describing for example a new data set, you can put it in ``lib/context`` 
to make it available to other people. In this case, make sure to update the list of available 
contexts :ref:`sec_context_list`.


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

 * Indent by tabs (make sure your editor does not sliently convert them to spaces). As
   indention is part of the python syntax it's important not to mix indention patters!
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


