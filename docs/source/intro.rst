.. _intro:

Introduction
============

About this project
------------------

Proteus is a scalable and versatile fingerprint and machine-learning code is a scalable and versatile machine-learning code. It collects several different algorithms that can be used to create fingerprints, perform dimensionality reduction and fit, atomistic and finite element calculations.

Proteus as of now is thought as standalone code. However, we aim to provide enough flexibility to interface it with other codes such as LAMMPS and PLUMED-2.0. It can be used as a C++ library as well as a python module. To be able to call it from python, we have used the pybind11 library.

Although at the moment is a serial-only code, we aim to write it in MPI so that it will be possible to take advantage of parallelization to speed up the calculations significantly.

It comes with a GNU lesser general public license of version 3, which means that it can be modified and freely distribute, although we take no responsibility for its misuse.


Compile the Code
----------------

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

.. _code_structure:

Structure of the Code
---------------------

The code is divided mainly in two parts: a pure C++ part and a python-binding interface.

The subroutines of the code that performs the expensive part of the calculation are written in C++14, and are collected in the **/src/** directory, as in the case of the :cpp:func:`cdist <cdist>`. This folder is completely agnostic of the python binding, and it should be kept in this way.

The python-bindings is obtained through Pybind11, and the binding subroutines are included in the **/bindings/** folder. Here, the bind_py_module.cc is the file that contains the main binding of the C++ package (in other words, is what is needed to use the syntax ``import Proteus``). The binding for each submodule of Proteus and its members are collected in files named bind_py_METHOD.cc. We decided to employ Pybind11 because of its seamless integration between Eigen and numpy. This allows the developer (and the user) to code fast and efficient algorithms easily, without losing the power of the C++ linear algebra as well as numpy simplicity. For more reference, please consult `Eigen interface from pybind11 <http://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html?highlight=eigen#pass-by-reference>`_.
