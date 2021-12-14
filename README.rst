librascal
================================================================================

.. start-intro

librascal is a versatile and scalable fingerprint and machine learning code. It
focuses on the efficient construction of representations of atomic structures,
that can then be fed to any supervised or unsupervised learning algorithm.
Simple regression code will be included for testing purposes, but the long-term
goal is to develop a separate collection of tools to this end.

librascal is currently considered a standalone code. However, we aim to provide
enough flexibility to interface it with other codes such as LAMMPS and
PLUMED-2.0. It can be used as a C++ library as well as a python module. To be
able to call it from python, we have used the pybind11 library.

Although at the moment is a serial-only code, we aim to write it in MPI so that
it will be possible to take advantage of parallelization to speed up the
calculations significantly. Parallelization is possible especially over atoms in
a structure (for large structures), over structures in a collection (for large
collections of small structures), or over components of a representation (for
representations with a large number of independent functions or components).

It comes with a GNU Lesser General Public License of version 3, which means that
it can be modified and freely distributed, although we take no responsibility
for its misuse.

For more information, have a look at the documentation_!

.. _documentation: https://cosmo-epfl.github.io/librascal/

Development
--------------------------------------------------------------------------------

The code is currently in the beta development phase, therefore we cannot
guarantee that the interface and data formats will not change, but it has been
in use for at least a year. Feedback and bug reports are welcome, as long as you
keep the above in mind.

.. end-intro

See `Helpers for Developers`_ below for some essential tools if you want to help
develop librascal. Be sure to also read `CONTRIBUTING.rst
<https://github.com/cosmo-epfl/librascal/blob/master/CONTRIBUTING.rst>`_ if you
plan on making a contribution.

Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before installing librascal, please make sure you have at least the following
packages installed:

+-------------+--------------------+
| Package     | Required version   |
+=============+====================+
| cmake       | 3.1 or higher      |
+-------------+--------------------+
| python      | 3.6 or higher      |
+-------------+--------------------+

To compile the code it is necessary to have a C++ compiler supporting C++14. For
gcc and clang this corresponds to versions

+-------------+--------------------+
| Compiler    | Required version   |
+=============+====================+
| gcc (g++)   | 4.9 or higher      |
+-------------+--------------------+
| clang       | 4.0 or higher      |
+-------------+--------------------+

Other necessary packages Eigen, PyBind11 and wigxjpf are downloaded
automatically when compiling librascal.


Installation
--------------------------------------------------------------------------------

.. start-install

When the dependencies are met, the python package can be installed with:

.. code:: bash

   pip install librascal

For optional features we have python packages required for the optional
librascal features:

+--------------------------+-------------+--------------------+
| Librascal feature        | Package     | Required version   |
+==========================+=============+====================+
| Feature compression      | skcosmo     | 0.1.0 or later     |
+--------------------------+-------------+--------------------+
| Rotational algebra       | sympy       | 1.4 or later       |
| (Clebsch-Gordan coeffs.) |             |                    |
+--------------------------+-------------+--------------------+

The dependencies for the optional features and the `introductive examples
<https://cosmo-epfl.github.io/librascal/examples/examples.html>`_ can be
installed with

.. code:: bash

   pip install -r requirements/common.txt

Compiling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To configure and compile the code with the default options, on \*nix systems
(Windows is not supported):

.. code:: shell

   mkdir build
   cd build
   cmake ..
   make

Customizing the build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The library supports several alternative builds that have additional
dependencies. Note that the ``ncurses`` GUI for cmake (ccmake) is quite helpful
to customize the build options.

Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Librascal source code is extensively tested (both c++ and python). The BOOST
unit_test_framework is required to build the tests (see BOOST.md for further
details on how to install the boost library). To build and run the tests:

.. code:: shell

   cd build
   cmake -DBUILD_TESTS=ON ..
   make
   ctest -V

You can also run the tests with Valgrind (a memory-error checker) by passing
``-DRASCAL_TESTS_USE_VALGRIND=ON`` to ``cmake``.

In addition to testing the behaviour of the code, the test suite also check for
formatting compliance with clang-format 8.0 or higher and black packages (these
dependencies are optional). To install these dependencies on Ubuntu:

.. code:: shell

   sudo apt-get install clang-format-8
   pip install -r requirements/testing.txt

Build Type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several build types are available Release (default), Debug and RelWithDebInfo.
To build an alternative mode

.. code:: shell

   cd build
   cmake -DCMAKE_BUILD_TYPE=Debug
   ..
   make

Or

.. code:: shell

   cd build
   cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo  \\
      CMAKE_C_FLAGS_RELWITHDEBUBINFO="-03 -g -DNDEBUG" ..
   make

Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The documentation relies on the sphinx (with nbsphinx and breathe extensions),
doxygen, pandoc, and graphviz packages. To install them on ubuntu:

.. code:: shell

  pip install -r requirements/doc.txt
  sudo apt-get install pandoc doxygen graphviz

Then to build the documentation run:

.. code:: shell

  cd build
  cmake -DBUILD_DOC=ON ..
  make doc

and open `build/docs/html/index.html` in a browser.

Bindings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Librascal relies on the pybind11 library to automate the generation of the
python bindings which are built by default. Nevertheless, to build only the c++
library:

.. code:: shell

   cd build
   cmake -DBUILD_BINDINGS=OFF ..
   make
