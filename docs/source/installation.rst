.. _installation

Installation
============


Prequisties
###########

Bofore installing Proteus, please make sure you have at least the following packages installed:

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

Other necessary packages (such as Eigen and PyBind11) are downloaded automatically whrn compiling Proteus.


Get Proteus
###########

It is strongly recommended to clone Proteus from its repository in GitHub.

.. code-block:: bash

    git clone LINK TO PROTEUS REPOSITORY

Build Proteus
#############

To configure and compile the code on UNIX, follow the command list:

.. code-block:: bash

  mkdir build
  cd build
  cmake ..
  make

Proteus checks  automatically for the optimal configuration corresponding to your system.


To Compile the documentation from the build folder:

.. code-block:: bash

  make dev_doc

Sphinx and breathe packages are required for building the documentation.


Run Proteus
###########

In order to run Proteus, you need to import the library into a Python code:

.. code-block:: python
    
    import Proteus as P


Advanced options
################

It is possible to link Proteus with other scientific calculation packages, s.a LAMMPS, ASE,... etc. A specific flag needs to be specified when building Proteus.
