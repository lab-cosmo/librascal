.. _installation

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
| python      |  3.4 or higher     |
+-------------+--------------------+

Other necessary packages (such as Eigen and PyBind11) are downloaded automatically whrn compiling Rascal.


Get Rascal
###########

It is strongly recommended to clone Rascal from its repository in GitHub.

.. code-block:: bash

    git clone LINK TO RASCAL REPOSITORY

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


To Compile the documentation from the build folder:

.. code-block:: console

  make dev_doc

Since Eigen depends on Mercurial, it may fail to download if you don't have the necessary dependencies. In that case, it may be sufficient to fix the dependencies or proceed by yourself.


Run Rascal
###########

In order to run Rascal, you need to import the library into a Python code:

.. code-block:: python
    
    import Rascal as P


Advanced options
################

It is possible to link Rascal with other scientific calculation packages, s.a LAMMPS, ASE,... etc. A specific flag needs to be specified when building Rascal.
