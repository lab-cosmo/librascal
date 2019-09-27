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
LAMMPS and PLUMED-2.0.  It can be used as a C++ library as well as a
python module.  To be able to call it from python, we have used the
pybind11 library.

Although at the moment is a serial-only code, we aim to write it in MPI
so that it will be possible to take advantage of parallelization to
speed up the calculations significantly.  Parallelization is possible especially
over atoms in a structure (for large structures), over structures in a
collection (for large collections of small structures), or over components of a
representation (for representations with a large number of independent functions
or components).

It comes with a GNU Lesser General Public License of version 3, which
means that it can be modified and freely distributed, although we take
no responsibility for its misuse.

Development
-----------

The code is currently in the alpha development phase; it is not yet
suitable for public use. Nevertheless, there is a significant amount of
functionality (including two tutorials) currently working and available
to test if you’re feeling adventurous. Feedback and bug reports are
welcome, as long as you keep the above in mind.

.. end-intro

Installation
------------

.. start-install

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
| ASE         | 3.18 or higher     |
+-------------+--------------------+

Other necessary packages (such as Eigen and PyBind11) are downloaded
automatically when compiling Rascal.

The following packages are required for building some optional features:

+------------------+-------------+--------------------+
| Feature          | Package     | Required version   |
+==================+=============+====================+
| Documentation    | pandoc      | (latest)           |
+------------------+-------------+--------------------+
|                  | sphinx      | 2.1.2              |
+------------------+-------------+--------------------+
|                  | breathe     | 4.13.1             |
+------------------+-------------+--------------------+
|                  | nbsphinx    | (latest)           |
+------------------+-------------+--------------------+

Compiling
~~~~~~~~~

To compile the code it is necessary to have CMake 3.0 and a C++ compiler
supporting C++14. During the configuration, it will automatically try to
download the external libraries on which it depends:

-  Eigen
-  Pybind11
-  Boost (only the unit test framework library)
-  Python3

And the following libraries to build the documentation:

-  Doxygen
-  Sphinx
-  Breathe

Beware, Python3 is mandatory. The code won’t work with a Python version
older than 3.

Using the package manager of your choice this yaml script should install all
required python packages required for the usage and development of rascal.

.. code:: yaml

    name: librascal-dev-env
    dependencies:
      - python=3.6 
      - pip:
        - numpy
        - matplotlib
        - scipy
        - mpmath
        - ase
        - ubjson
        - cpplint
        - sphinx=2.1.2
        - sphinx_rtd_theme
        - breathe=4.13.1
        - nbsphinx
        - jupyter
        - qml
        - autopep8
        - pytest
        - tqdm

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

1. Tests

   Librascal source code is extensively tested (both c++ and python).
   The BOOST unit_test_framework is requiered to build the tests (see
   BOOST.md for further details on how to install the boost library). To
   build and run the tests:

   .. code:: shell

      cd build
      cmake -DBUILD_TESTS=ON ..
      make
      ctest -V

   In addition to testing the behaviour of the code, the test suite also check
   for formatting compliance with the clang-format and autopep8 packages (these
   dependencies are optional). To install these dependencies on ubuntu:

   .. code:: shell

      sudo apt-get install clang-format
      pip3 install autopep8

2. Build Type

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

3. Documentation

   The documentation relies on the sphinx (with nbsphinx and breathe
   extensions), doxygen, pandoc, and graphviz
   packages. To install them on ubuntu:

   .. code:: shell

     pip3 install sphinx sphinx_rtd_theme breathe nbsphinx
     sudo apt-get install pandoc doxygen graphviz

   Then to build the documentation run:

   .. code:: shell

     cd build
     cmake -DENABLE_DOC=ON ..
     make doc

   and open :file:`build/docs/html/intro.html` in a browser.  The
   Doxygen html build is still available under
   :file:`build/docs/dox_html/index.html` but it is being phased out in favour
   of Sphinx/breathe.

4. Helpers for Developers

   -  To remove all the cmake files/folders except for the external
      library (enable glob and remove):

   .. code:: shell

      shopt -s extglob
      rm -fr -- !(external|third-party)

   -  To help developers conform their contribution to the coding
      convention, the formatting of new functionalities can be automated
      using clang-format (for the c++ files) and autopep8 (for the
      python files). The .clang-format and .pycodestyle files define
      common settings to be used.

      To enable these functionalities (optional) you can install these
      tools with:

      .. code:: shell

         sudo apt-get install clang-format
         pip install autopep8

      The automatic formating of the c++ and python files can be
      trigered by:

      .. code:: shell

         cd build
         cmake ..
         make pretty-cpp
         make pretty-python

      Please use these tools with caution as they can potentially
      introduce unwanted changes to the code. If code needs to be
      specifically excluded from auto formatting, e.g. a matrix which
      should be human-readable, code comments tells the formatters to
      ignore lines:

      C++

      .. code:: C++

         // clang-format off
         SOME CODE TO IGNORE
         // clang-format on

      python

      .. code:: python

         SOME LINE TO IGNORE # noqa

      where ``noqa`` stands for ``no`` ``q``\ uality ``a``\ ssurance.

5. Bindings

   Librascal relies on the pybind11 library to automate the generation
   of the python bindings which are built by default. Nevertheless, to
   build only the c++ library:

   .. code:: shell

      cd build
      cmake -DBUILD_BINDINGS=OFF ..
      make

6. Installing rascal

   .. code:: shell

      mkdir ../build
      cd build
      cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DBUILD_BINDINGS=ON ..
      make install

Miscellaneous Information
-------------------------

-  Common cmake flags:

   -  -DCMAKE_C_COMPILER
   -  -DBUILD_BINDINGS
   -  -DUSER
   -  -DINSTALL_PATH
   -  -DCMAKE_BUILD_TYPE
   -  -DENABLE_DOC
   -  -DBUILD_TESTS

-  Special flags:

   -  -DBUILD_BINDINGS:

      -  ON (default) -> build python binding
      -  OFF -> does not build python binding

   -  -DINSTALL_PATH:

      -  empty (default) -> does not install in a custom folder
      -  custom string -> root path for the installation

   -  -DUSER:

      -  OFF (default) -> changes nothing
      -  ON -> install root is in the user’s home directory, i.e.
         ``~/.local/``

To build librascal as a docker environment:

.. code:: shell

   sudo docker build -t test -f ./docker/install_env.dockerfile  .
   sudo docker run -it -v /path/to/repo/:/home/user/  test
