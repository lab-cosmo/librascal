.. Rascal documentation master file, created by
   sphinx-quickstart on Thu Mar  1 13:58:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Rascal's documentation!
===================================

Rascal is a scalable and versatile fingerprint and machine-learning code is a scalable and versatile machine-learning code. It collects several different algorithms that can be used to create fingerprints, perform dimensionality reduction and fit, atomistic and finite element calculations.

Rascal as of now is thought as standalone code. However, we aim to provide enough flexibility to interface it with other codes such as LAMMPS and PLUMED-2.0. It can be used as a C++ library as well as a python module. To be able to call it from python, we have used the pybind11 library.

Although at the moment is a serial-only code, we aim to write it in MPI so that it will be possible to take advantage of parallelization to speed up the calculations significantly.

It comes with a GNU Lesser General Public License of version 3, which means that it can be modified and freely distribute, although we take no responsibility for its misuse.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   intro
   whitepaper
   installation
   tutorials
   representations
   models
   neighbor-list-manager
   developer
   auto

Indices and tables
==================

* :ref:`intro`
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



