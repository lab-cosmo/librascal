.. Rascal documentation master file, created by
   sphinx-quickstart on Thu Mar  1 13:58:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to librascal's documentation!
=====================================

.. include:: ../../README.rst
   :start-after: start-intro
   :end-before: end-intro


.. role:: python(code)
   :language: python

.. _code_structure:

Structure of the Code
---------------------

The code is divided mainly in two parts: a pure C++ part and a python-binding interface.

The subroutines of the code that performs the expensive part of the calculation
are written in C++14, and are collected in the :file:`/src/` directory, as in the
case of the :cpp:func:`cdist <cdist>`. This folder is completely agnostic of the
python bindings, and it should be kept in this way.

The python-bindings is obtained through Pybind11, and the binding subroutines
are included in the :file:`/bindings/` folder. Here, the bind_py_module.cc is
the file that contains the main binding of the C++ package (in other words, is
what is needed to use the syntax :python:`import librascal`). The binding for
each submodule of Rascal and its members are collected in files named
bind_py_METHOD.cc. We decided to employ Pybind11 because of its seamless
integration between Eigen and numpy. This allows the developer (and the user) to
code fast and efficient algorithms easily, without losing the power of the C++
linear algebra as well as numpy simplicity. For more reference, please consult
`Eigen interface from pybind11 (outbound)
<http://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html?highlight=eigen#pass-by-reference>`_.

Developers
----------

Rascal is jointly developed between the `COSMO (outbound) <https://cosmo.epfl.ch>`_ and
`LAMMM (outbound) <https://lammm.epfl.ch>`_ groups at EPFL.  The active developers are
Michele Ceriotti, Federico Giberti, Alexander Gocsinski, Till Junge, Félix
Musil, Markus Stricker, Max Veit, and Michael Willatt.

This documentation has been contributed by Chiheb Ben Mahmoud, Michele Ceriotti,
Federico Giberti, Klim Goldshtein, Till Junge, Markus Stricker, Félix Musil,
and is currently maintained by Max Veit.


Contents
----------

.. toctree::
   :maxdepth: 2

   installation
   tutorials/index
   whitepaper
   SOAP
   representation/representations
   models/models
   io/io
   neighbour_list/neighbour-list-manager
   dev_guide/developer
   reference/index
   bibliography


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
