librascal
=========

.. start-intro

librascal is a versatile and scalable fingerprint and machine learning
code. It focuses on the efficient construction of representations of
atomic structures, that can then be fed to any supervised or
unsupervised learning algorithm. Simple regression code will be included
for testing purposes, but the long-term goal is to develop a separate
collection of tools to this end.

librascal is currently considered a standalone code. However, we aim to
provide enough flexibility to interface it with other codes such as
LAMMPS and PLUMED-2.0. It can be used as a C++ library as well as a
python module. To be able to call it from python, we have used the
pybind11_ library.

Although at the moment is a serial-only code, we aim to write it in MPI
so that it will be possible to take advantage of parallelization to
speed up the calculations significantly. Parallelization is possible especially
over atoms in a structure (for large structures), over structures in a
collection (for large collections of small structures), or over components of a
representation (for representations with a large number of independent functions
or components).

librascal is distributed under the GNU Lesser General Public License version 2.1
or later, at your convenience. This means that it can be modified and freely
distributed, although we take no responsibility for its misuse.

For more information, have a look at the documentation_!

.. _documentation: https://lab-cosmo.github.io/librascal/
.. _pybind11: https://pybind11.readthedocs.io

Development
-----------

The code is currently in the alpha development phase; it is not yet
suitable for public use. Nevertheless, there is a significant amount of
functionality (including two tutorials) currently working and available
to test if you’re feeling adventurous. Feedback and bug reports are
welcome, as long as you keep the above in mind.

.. end-intro

See `Helpers for Developers`_ below for some essential tools if you want to help
develop libRascal.  Be sure to also read `CONTRIBUTING.rst <CONTRIBUTING.rst>`_
if you plan on making a contribution.

Installation
------------

.. start-install

The installation of the library for python use can be done simply with:

.. code:: bash

   pip install .

assuming that `python` 3.5 (or higher) and `gcc` or `clang` are available.

Dependencies
~~~~~~~~~~~~

Before installing librascal, please make sure you have at least the
following packages installed:

+-------------+--------------------+
| Package     | Required version   |
+=============+====================+
| gcc (g++)   | 4.9 or higher      |
+-------------+--------------------+
| clang       | 4.0 or higher      |
+-------------+--------------------+
| cmake       | 2.8 or higher      |
+-------------+--------------------+
| python      | 3.6 or higher      |
+-------------+--------------------+
| numpy       | 1.13 or higher     |
+-------------+--------------------+
| scipy       | 1.4.0 or higher    |
+-------------+--------------------+
| ASE         | 3.18 or higher     |
+-------------+--------------------+

Other necessary packages (such as Eigen and pybind11) are downloaded
automatically when compiling Rascal.

The following packages are required for some optional features:

+--------------------------+-------------+--------------------+
| Feature                  | Package     | Required version   |
+==========================+=============+====================+
| Feature compression      | skcosmo     | 0.1.0 or later     |
+--------------------------+-------------+--------------------+
| Rotational algebra       | sympy       | 1.4 or later       |
| (Clebsch-Gordan coeffs.) |             |                    |
+--------------------------+-------------+--------------------+
| Building documentation   | pandoc      | (latest)           |
+--------------------------+-------------+--------------------+
|                          | sphinx      | 2.1.2 or later     |
+--------------------------+-------------+--------------------+
|                          | breathe     | 4.14.1 or later    |
+--------------------------+-------------+--------------------+
|                          | nbsphinx    | 0.8.1 or later     |
+--------------------------+-------------+--------------------+

Compiling
~~~~~~~~~

To compile the code it is necessary to have CMake 3.0 and a C++ compiler
supporting C++14. During the configuration, it will automatically try to
download the external libraries on which it depends:

-  Eigen
-  pybind11
-  Boost (only the unit test framework library)
-  Python3

And the following libraries to build the documentation:

-  Doxygen
-  Sphinx
-  Breathe

Beware, Python3 is mandatory. The code won’t work with a Python version
older than 3.

You can then use pip to install all python packages required for the usage
and development of rascal:

.. code:: bash

    pip install -r requirements.txt

To configure and compile the code with the default options, on \*nix
systems (Windows is not supported):

.. code:: shell

   mkdir build
   cd build
   cmake ..
   make

Customizing the build
~~~~~~~~~~~~~~~~~~~~~

The library supports several alternative builds that have additional
dependencies. Note that the ``ncurses`` GUI for cmake (ccmake) is quite
helpful to customize the build options.

Tests
^^^^^

Librascal source code is extensively tested (both c++ and python).
The BOOST unit_test_framework is required to build the tests (see
BOOST.md for further details on how to install the boost library). To
build and run the tests:

.. code:: shell

   cd build
   cmake -DBUILD_TESTS=ON ..
   make
   ctest -V

You can also run the tests with Valgrind (a memory-error checker) by passing
``-DRASCAL_TESTS_USE_VALGRIND=ON`` to ``cmake``.

In addition to testing the behaviour of the code, the test suite also check
for formatting compliance with clang-format 8.0 or higher and black packages
(these dependencies are optional). To install these dependencies on Ubuntu:

.. code:: shell

   sudo apt-get install clang-format-8
   pip3 install black

Build Type
^^^^^^^^^^

Several build types are available Release (default), Debug and
RelWithDebInfo. To build an alternative mode

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
^^^^^^^^^^^^^

The documentation relies on the sphinx (with nbsphinx and breathe
extensions), doxygen, pandoc, and graphviz
packages. To install them on ubuntu:

.. code:: shell

  pip3 install sphinx sphinx_rtd_theme breathe nbsphinx
  sudo apt-get install pandoc doxygen graphviz

Then to build the documentation run:

.. code:: shell

  cd build
  cmake -DBUILD_DOC=ON ..
  make doc

and open :file:`build/docs/html/index.html` in a browser.

Bindings
^^^^^^^^

Librascal relies on the pybind11 library to automate the generation
of the python bindings which are built by default. Nevertheless, to
build only the c++ library:

.. code:: shell

   cd build
   cmake -DBUILD_BINDINGS=OFF ..
   make

Installing rascal
^^^^^^^^^^^^^^^^^
To install the python library with c++ bindings:

.. code:: shell

   pip install .


Helpers for Developers
~~~~~~~~~~~~~~~~~~~~~~

Deepclean
^^^^^^^^^

To remove all the cmake files/folders except for the external
library (enable glob and remove):

.. code:: shell

   shopt -s extglob
   rm -fr -- !(external|third-party)

Automatic code formatting
^^^^^^^^^^^^^^^^^^^^^^^^^

To help developers conform their contribution to the coding
convention, the formatting of new functionalities can be automated
using clang-format (for the c++ files) and black (for the
python files). The .clang-format and .pycodestyle files define
common settings to be used.

To enable these functionalities (optional) you can install these
tools with:

.. code:: shell

   sudo apt-get install clang-format
   pip install black

The automatic formatting of the c++ and python files can be
triggered by:

.. code:: shell

   cd build
   cmake ..
   make pretty-cpp
   make pretty-python

Please use these tools with caution as they can potentially
introduce unwanted changes to the code. If code needs to be
specifically excluded from auto formatting, e.g. a matrix which
should be human-readable, code comments tells the formatters to
ignore lines:

- C++

  .. code:: C++

     // clang-format off
     SOME CODE TO IGNORE
     // clang-format on

- python

  .. code:: python

     SOME LINE TO IGNORE # noqa

  where ``noqa`` stands for ``no`` ``q``\ uality ``a``\ ssurance.

Jupyter notebooks
^^^^^^^^^^^^^^^^^

If you are contributing any code in IPython/Jupyter notebooks, *please*
install the `nbstripout` extension (available e.g. on
`github <https://github.com/kynan/nbstripout#installation>`_ and
`PyPI <https://pypi.org/project/nbstripout/>`_).  After installing,
activate it for this project by running:

.. code:: shell

   nbstripout --install --attributes .gitattributes

from the top-level repository directory.  Please note that that
``nbstripout`` will not strip output from cells with the metadata fields
``keep_output`` or ``init_cell`` set to ``True``, so use these fields
judiciously.  You can ignore these settings with the following command:

.. code:: shell

   git config filter.nbstripout.extrakeys '\
      cell.metadata.keep_output cell.metadata.init_cell'

(The keys ``metadata.kernel_spec.name`` and
``metadata.kernel_spec.display_name`` may also be useful to reduce diff
noise.)

Nonetheless, it is highly discouraged to contribute code in the form of
notebooks; even with filters like ``nbstripout`` they're a hassle to use
in version control.  Use them only for tutorials or *stable* examples that
are either meant to be run *interactively* or are meant to be processed by
`sphinx` (`nbsphinx <https://nbsphinx.readthedocs.io/en/latest/>`_) for
inclusion in the
`tutorials page <https://lab-cosmo.github.io/librascal/tutorials/tutorials.html>`_.

Miscellaneous Information
-------------------------

-  Common cmake flags:

   -  -DCMAKE_CXX_COMPILER
   -  -DCMAKE_C_COMPILER
   -  -DCMAKE_BUILD_TYPE
   -  -DBUILD_BINDINGS
   -  -DINSTALL_PATH
   -  -DBUILD_DOC
   -  -DBUILD_TESTS

-  Special flags:

   -  -DBUILD_BINDINGS:

      -  ON (default) -> build python binding
      -  OFF -> does not build python binding

   -  -DINSTALL_PATH:

      -  empty (default) -> does not install in a custom folder
      -  custom string -> root path for the installation


To build librascal as a docker environment:

.. code:: shell

   sudo docker build -t test -f ./docker/install_env.dockerfile  .
   sudo docker run -it -v /path/to/repo/:/home/user/  test
