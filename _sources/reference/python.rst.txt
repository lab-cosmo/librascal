.. _auto_python:


Python
------

This list is incomplete. You can help by *expanding it*!

Atoms and cells
===============

Some base functionality for working with structures and collections thereof in
``librascal``

.. autoclass:: rascal.neighbourlist.structure_manager.AtomsList

.. autofunction:: rascal.neighbourlist.structure_manager.mask_center_atoms_by_id

.. autofunction:: rascal.neighbourlist.structure_manager.mask_center_atoms_by_species

Representations
===============

Representations are the primary classes in ``librascal`` used to compute
structural representations (features, descriptors) from a list of atoms.

.. autoclass:: rascal.representations.SortedCoulombMatrix
   :members:

.. autoclass:: rascal.representations.SphericalInvariants
   :members:

.. autoclass:: rascal.representations.SphericalCovariants
   :members:

Models
======

Also included is an optimized implementation of kernel ridge regression (KRR,
also equivalent to Gaussian approximation potentials aka GAP).  Both fitting and
evaluating of models is implemented; the evaluation in particular is optimized
for use in MD simulations.

.. autoclass:: rascal.models.Kernel
   :members:

.. autoclass:: rascal.models.KRR
   :members:

.. autofunction:: rascal.models.compute_KNM


.. _filters:

Filters
^^^^^^^

Filters are used mainly to select rows and columns from feature (or kernel)
matrices in order to reduce their dimensionality and make the fitting problem
tractable, or just more efficient.  The following filters, implemented in
`scikit-cosmo (outbound) <https://scikit-cosmo.readthedocs.io/en/latest/selection.html>`_, are available:

.. autoclass:: rascal.utils.FPSFilter
    :members:

.. autoclass:: rascal.utils.CURFilter
    :members:

Both inherit the interface of the following base class:

.. autoclass:: rascal.utils.filter.Filter


IO
===

Utilities for loading and saving Rascal objects (especially models)

.. autofunction:: rascal.utils.dump_obj

.. autofunction:: rascal.utils.load_obj

.. autoclass:: rascal.utils.BaseIO
   :members:


