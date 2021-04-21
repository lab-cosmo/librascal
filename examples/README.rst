Examples of using `librascal`
=============================

This folder contains hands-on examples and tutorials that demonstrate how to use librascal for atomistic machine learning.  They currently focus on the task of fitting a model and using it to run machine learning molecular dynamics simulations (ML-MD).  The two notebooks that currently demonstrate this are:

- :file:`MLIP_example.ipynb`: Generic example of fitting a machine learning potential and running it using `ASE <https://wiki.fysik.dtu.dk/ase>`_
- :file:`zundel_i-PI.ipynb` (with data in the :file:`i-PI/zundel/` directory): Example of fitting a potential and running using the `i-PI <http://ipi-code.org/>`_ code

These notebooks are small, fast, and lightweight, and are automatically checked to make sure they run with the current version of the code.

Further examples can be found in the various directories, though they are not guaranteed to be working or up-to-date.  They are not recommended for new users learning to use the library.

- :file:`python/`: Examples of using the Python bindings, e.g. for model training
- :file:`cpp/`: Examples of using the low-level C++ interface, can be compiled by adding ``-DBUILD_EXAMPLES`` to the CMake command line
- :file:`needs_updating/`: Notebooks that need updating to work with the most recent version of the code; they may still contain snippets useful to experienced users
