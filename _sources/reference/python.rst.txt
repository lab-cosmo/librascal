.. _auto_python:


Python
------

This list is incomplete. You can help by *expanding it*!

Atoms and cells
===============

.. autoclass:: rascal.neighbourlist.structure_manager.AtomsList

.. autofunction:: rascal.neighbourlist.structure_manager.mask_center_atoms_by_id

.. autofunction:: rascal.neighbourlist.structure_manager.mask_center_atoms_by_species

Representations
===============

.. autoclass:: rascal.representations.SortedCoulombMatrix
   :members:

.. autoclass:: rascal.representations.SphericalInvariants
   :members:

.. autoclass:: rascal.representations.SphericalCovariants
   :members:

Models
======
.. autoclass:: rascal.models.Kernel
   :members:

.. autoclass:: rascal.models.KRR
   :members:

.. autofunction:: rascal.models.compute_KNM

IO
===

.. autofunction:: rascal.utils.dump_obj

.. autofunction:: rascal.utils.load_obj

.. autoclass:: rascal.utils.BaseIO
   :members:


