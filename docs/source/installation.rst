.. _installation:

Installation
============


Prequisties
###########

Bofore installing Rascal, please make sure you have at least the following packages installed:

+-------------+--------------------+
| Package     |  Required version  |
+=============+====================+
| openmpi     |  1.8 or higher     |
+-------------+--------------------+
| gcc or Clang|  4.9 or higher     |
+-------------+--------------------+
| cmake       |  2.8 or higher     |
+-------------+--------------------+
| python      |  3.6 or higher     |
+-------------+--------------------+
| numpy       |  1.13 or higher    |
+-------------+--------------------+


Other necessary packages (such as Eigen and PyBind11) are downloaded automatically when compiling Rascal.


Get Rascal
###########

It is strongly recommended to clone Rascal from its repository in GitHub.

.. code-block:: bash

    git clone https://github.com/cosmo-epfl/librascal.git

Build Rascal
#############

To compile the code it is necessary to have CMake 3.0 and a C++ compiler supporting C++14. During the configuration, it will automatically try to download the external libraries on which it depends:


- Eigen
- Pybind11
- Boost (only the unit test framework  library)
- Doxygen
- Sphinx
- Breathe
- Python3

Beware, Python3 is mandatory. The code won't work with a Python version older than 3.


To configure and compile the code on UNIX, follow the command list:

.. code-block:: console

  mkdir build
  cd build
  cmake ..
  make


Alternative build
#################

The library supports several alternative build that have additional depencies.
Note that the curse gui for cmake (ccmake) is quite helpfull to customize the build options.

1. Tests

Librascal source code is extensively tested (both c++ and python). The BOOST unit_test_framework is requiered to build the tests (see below for further details on how to install the boost library).
To build and run the tests:

.. code-block:: bash

    cd build
    cmake -DBUILD_TESTS=ON ..
    make
    ctest -V


In addition to testing the behaviour of the code, the test suite also check for formatting compliance with the clang-format and autopep8 packages (these dependencies are optional).
To install these dependencies on ubuntu:

.. code-block:: bash

    sudo apt-get install clang-format
    pip3 install autopep8


2. Build Type

Several build types are available Release (default), Debug and RelWithDebInfo. To build an alternative mode:

.. code-block:: bash

    cd build
    cmake -DCMAKE_BUILD_TYPE=Debug  ..
    make


Or

.. code-block:: bash

    cd build
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo  CMAKE_C_FLAGS_RELWITHDEBUBINFO="-03 -g -DNDEBUG" ..
    make


3. Documentation

The current documentation relies on the doxygen package. To install it on ubuntu:

.. code-block:: bash

    sudo apt-get install doxygen


The documentation is located in the librascal/docs/documentation/html folder. The source files for the documentation are located in the librascal/docs/src folder. 
To rebuild the documentation, run the 

.. code-block:: bash

    make dev_doc

in the build folder.

4. Bindings

Librascal relies on the pybind11 library to automate the generation of the python bindings which are built by default. Nevertheless, to build only the c++ library:

.. code-block:: bash

    cd build
    cmake -DBUILD_BINDINGS=OFF ..
    make


5. Helpers for Developers

* To remove all the cmake files/folders except for the external library (enable glob and remove):

.. code-block:: bash

    shopt -s extglob
    rm -fr -- !(external|third-party)


* To help developers conform their contribution to the coding convention, the formating of new functionalities can be automated using clang-format (for the c++ files) and autopep8 (for the python files). The .clang-format and .pycodestyle files define common settings to be used.

* To enable these functionalities (optional) you can install these tools with:

.. code-block:: bash

    sudo apt-get install clang-format
    pip install autopep8


The automatic formating of the c++ and python files can be trigered by:

.. code-block:: bash

    cd build
    cmake ..
    make pretty-cpp
    make pretty-python


Please use these tools with caution as they can potentially introduce unwanted changes to the code.
If code needs to be specifically excluded from auto formatting, e.g. a matrix which should be human-readable, code comments tells the formatters to ignore lines:

* C++

.. code-block:: bash

    // clang-format off
    SOME CODE TO IGNORE
    // clang-format on
    

* Python

.. code-block:: bash

    SOME LINE TO IGNORE # noqa


where <b>`noqa`</b> stands for <b>no</b> <b>q</b>uality <b>a</b>ssurance.

Misceleaneous Information
########################

* Common cmake flag:
  + -DCMAKE_C_COMPILER
  + -DBUILD_BINDINGS
  + -DUSER
  + -DINSTALL_PATH
  + -DCMAKE_BUILD_TYPE
  + -DENABLE_DOC
  + -DBUILD_TESTS

* Special flags:
  + -DBUILD_BINDINGS:
    + ON (default) -> build python binding
    + OFF -> does not build python binding
  + -DINSTALL_PATH:
    + empty (default) -> does not install in a custom folder
    + custom string -> root path for the installation
  + -DUSER:
    + OFF (default) -> changes nothing
    + ON -> install root is in the user's home directory, i.e. ~/.local/


To build libRascal with as docker environement:

.. code-block:: bash

    sudo docker build -t test -f ./docker/install_env.dockerfile  .
    sudo docker run -it -v /path/to/repo/:/home/user/  test



Run Rascal
###########

In order to run Rascal, you need to import the library into a Python code:

.. code-block:: python
    
    import Rascal as R


Advanced options
################

It is possible to link Rascal with other scientific calculation packages, s.a LAMMPS, ASE,... etc. A specific flag needs to be specified when building Rascal.
