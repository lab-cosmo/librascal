.. _neighbour-list-manager:

Neighbour List Manager
=====================

Librascal provides a user friendly, flexible and efficient infrastructure to build neighbour list meeting the specific needs of many represenations of atomic neighbourhood.
A neighbour list in librascal manages and/or builds all the data relevant to the atomic structure it refers to, e.g. the atomic structure (positions, atomic types, cell, periodic boundary conditions), the distances between a central atom and its neighbours, the representations associated to the atomic structure...
The main purpose of gathering all these informations in one place is to reduce as much as possible the access to these quantities with explicit indices (see :ref:`nl-for-user` for a detailed description of how to interact with them).

The purpose of neighbour lists in librascal is to provide the iteration patterns dictated by the definition of a representation from an input structure which could be a raw atomic structure or an already existing neighbour list provided by an external code such as LAMMPS. Clearly both the starting and ending state of a neighbour list differ greatly between representations and framework. To adress this issue librascal defines two family of objects: the :class:`StructureManagerX` and the :class:`AdaptorY` that follow the interface defined the :cpp:class:`StructureManager <rascal::StructureManager>` class (`X` and `Y` are names referring to the function of the object).
A structure manager is meant to handle the input structure like in :cpp:class:`StructureManagerCenters <rascal::StructureManagerCenters>` or :cpp:class:`StructureManagerLammps <rascal::StructureManagerLammps>` while an adaptor modifies the neighbour list so that it matches the requirements of a representation.



.. toctree::
    :caption: General descriptions
    :maxdepth: 2

    adaptors
    species-manager
