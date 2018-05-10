.. _intro:

Introduction
============

.. contents::
   :local:
   

.. _code_structure:

Structure of the Code
---------------------

The code is divided mainly in two parts: a pure C++ part and a python-binding interface.

The subroutines of the code that performs the expensive part of the calculation are written in C++14, and are collected in the **/src/** directory, as in the case of the :cpp:func:`cdist <cdist>`. This folder is completely agnostic of the python binding, and it should be kept in this way.

The python-bindings is obtained through Pybind11, and the binding subroutines are included in the **/bindings/** folder. Here, the bind_py_module.cc is the file that contains the main binding of the C++ package (in other words, is what is needed to use the syntax ``import Rascal``). The binding for each submodule of Rascal and its members are collected in files named bind_py_METHOD.cc. We decided to employ Pybind11 because of its seamless integration between Eigen and numpy. This allows the developer (and the user) to code fast and efficient algorithms easily, without losing the power of the C++ linear algebra as well as numpy simplicity. For more reference, please consult `Eigen interface from pybind11 <http://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html?highlight=eigen#pass-by-reference>`_.
