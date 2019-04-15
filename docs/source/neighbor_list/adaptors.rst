.. _adaptors:

Adaptors
~~~~~~~~

.. contents::
   :local:

Overview
********
Adaptors are objects which are *stacked* on top of a :class:`~StructureManager` (or ``AdatorX`) to provide a new (extended or sorted) view on an atomic structure or neighbour list. The name *adaptor* is chosen because they are used to adapt and prepare an atomic structure for a given calculation based on e.g. pair relationships.
Each adaptor preserves the interface for interacting with an atomic structure.

This can for example be neighbourlist, filter by strict cutoffs, increase a half neighbourlist to a full neighbourlist.

A typical usage scenario is an :class:`~AdaptorStrict` on top of an :class:`~AdaptorNeighbourList` on top of a :class:`~StructureManager`. By adapting the initial atomic structure in this way, the final object is an iterable strict neighbourlist.

Each adaptor only stores the additional information it generates during construction. No atoms or neighbourlists are stored twice. If the requested information does not exist at this `Layer`, it is piped to the underlying adaptor/manager.


AdaptorNeighbourList
********************
Since the purpose of this library is the efficient calculation of atomic representations, this adaptor is probably the most fundamental. Without a neighbourlist, there is no relationship of an atom between its environment.

The :class:`~AdaptorNeighbourList` provides the functionality to build a neighbourlist from a list of atoms, their cell and periodicity information. By default a full neighbourlist is built based on the cell linked-lists and triclinicity is accounted for.

The basic idea is to anchor a mesh at the origin of the supplied cell. This overlaid mesh is then repeated until it is large enough to have at least one cutoff length in all direction. Periodicity is now ensured by adding ghost atoms through shifting all center atoms by the supplied lattice vectors.
All center and ghost atoms are then sorted into boxes. These boxes are indexed and the neighbourhood for each atom is built based on it's surrounding.
The resulting neighbourlist is full and not strict.

One peculiarity has to be mentiond. It is the flag ``consider_ghost_neighbours``.
The standard behaviour of the adaptor is to provide neighbours of the initial list of atoms. It does not provide a neighbourlist for the ghost atoms.
This corresponds to ``consider_ghost_neighbours=false``. And it is fine if only a pair list is needed for whatever comes afterwards.
But if triplets are constructed subsequently, the neighbours of ghost atoms are needed. Some atomic representations have three-body terms and they can not be constructed if the ghost atoms do not have neighbours.

AdaptorHalfList
***************
:class:`~AdaptorHalfList` provides to functionality to reduce a full neighbourlist to a half list. After construction the iteration over all possible pairs in an atomic structure yields each pair only once.
The choice for including or not including a pair is based on the atomic indices of the participating atoms. Only pairs where the first (i-atom) index is smaller than the second (j-atom) index are included.
Atomic indices are rascal-internal and the user does not have control over atomic indices. This ensures that the condition is a safe check.

AdaptorFullList
***************
:class:`~AdaptorFullList` provides the functionality to blow up a half neighbourlist to a full neighbourlist.
After construction the iteration over pairs yields each atom at least once as an i-atom.
By construction the logic of layering in rascal for pairs is broken. Hence layer information for all pairs (also the ones that previously existed) is reset.

This means that a previously calculated :class:`~Property` is not accessible any more with :class:`~ClusterRef`s of this type. Or it is partly accessible, but meaningless since the layering is destroyed.

AdaptorStrict
*************
:class:`~AdaptorStrict` provides the functionality to only iterate over neighbours that are actually within a cutoff.
Can must be stacked on top of an existing neighbourlist.
During calculation of the strict neighbourlist the actual distance between neighbouring atoms has to be calculated.
The distance (often used for calculating representations) and direction (used in derivatives) is stored for later possible reuse.

Using this adaptor ensures that later calculation of (possibly very expensive) representations is limited to the atoms which are in each other's cutoff.

AdaptorMaxOrder
***************
:class:`~AdaptorMaxOrder` increases and existing neighbourlist by one Order. I.e. a pair list becomes a triplet list, a triplet list a quadruplet list, etc.
It must be stacked on an existing neighbourlist.

The functionality is provided by taking a ``Cluster`` of a given order and adding all neighbours (pairs) of the constituting atoms as neighbours of this cluster.

(Currently there is a peculiarity. Adding all neighbours results in an ambiguity. The triplet list is not strict any more. This is expected to be sorted out during the implementation of ``Representations``)

.. _`adaptor filter`:
AdaptorFilter
*************
:class:`~AdaptorFilter` is a pure virtual which can not be instantiated. An daughter class of it provides a filtered (masked) view on an existing :class:`~StructureManager` or :class:`~Adaptor`.

Any daughter class needs to implement the :meth:`~perform_filtering()` function. When used in this way, it could in principle substitute the :class:`~AdaptorStrict`.

Currently only the :class:`~SpeciesManager` uses this class. Within this context the :class:`~AdaptorFilter` is used as a bucket to store a species-sorted view. Please refer to its documentation for details :ref:`species-manager`.
