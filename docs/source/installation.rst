.. _installation:

Installation
============

Get librascal
-------------

It is strongly recommended to clone librascal from its repository in GitHub.

.. code-block:: bash

    git clone https://github.com/lab-cosmo/librascal.git


.. include:: ../../README.rst
   :start-after: start-install


Run Rascal
----------

In order to run Rascal, you need to import the library into a Python code:

.. code-block:: python

    import rascal
    from rascal.representations import *


Advanced options
----------------

It is possible to link Rascal with other scientific calculation packages, like
LAMMPS, ASE, i-PI, and n2p2.  These interfaces are still a work in progress.
