About this project
==================

A scalable and versatile fingerprint and machine-learning code.

Compile the Code
================

.. code-block:: console

  mkdir build
  cd build
  cmake ./..
  make

To Compile the documentation from the build folder:

.. code-block:: console

  make dev_doc


Structure of the Code
=====================

The essential and expensive procedures are done in c++14 and stored in the src directory like in :cpp:func:`cdist <cdist>`.
The pybind11 bindings are stored in the bindings folder. The bind_py_module.cc file collects all the binding generator functions defined in the like prot_dist_mat.cc .

For the ease of interfacing c++ and python, we use the `Eigen interface from pybind11 <http://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html?highlight=eigen#pass-by-reference>`_.

