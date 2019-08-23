Adaptors
==========

Overview
----------
Adaptors are objects which are stacked on top of a StructureManager (or ``AdaptorX`) to provide a new (extended or sorted) view on an atomic structure or neighbour list. The name "adaptor" is chosen because they are used to adapt and prepare an atomic structure for a given calculation based on e.g. pair relationships. Each adaptor preserves the interface for interacting with an atomic structure.

This can for example be neighbourlist, filter by strict cutoffs, increase a half neighbourlist to a full neighbourlist.

A typical usage scenario is an AdaptorStrict on top of an AdaptorNeighbourList on top of a StructureManager. By adapting the initial atomic structure in this way, the final object is an iterable strict neighbourlist.

Each adaptor only stores the additional information it generates during construction. No atoms or neighbourlists are stored twice. If the requested information does not exist at this Layer, it is piped to the underlying adaptor/manager.
AdaptorNeighbourList
--------------------
Since the purpose of this library is the efficient calculation of atomic representations, this adaptor is probably the most fundamental. Without a neighbourlist, there is no relationship of an atom between its environment.

The AdaptorNeighbourList provides the functionality to build a neighbourlist from a list of atoms, their cell and periodicity information. By default a full neighbourlist is built based on the cell linked-lists and triclinicity is accounted for.

The basic idea is to anchor a mesh at the origin of the supplied cell. This overlaid mesh is then repeated until it is large enough to have at least one cutoff length in all direction. Periodicity is now ensured by adding ghost atoms through shifting all center atoms by the supplied lattice vectors. All center and ghost atoms are then sorted into boxes. These boxes are indexed and the neighbourhood for each atom is built based on it’s surrounding. The resulting neighbourlist is full and not strict.

One peculiarity has to be mentiond. It is the flag consider_ghost_neighbours. The standard behaviour of the adaptor is to provide neighbours of the initial list of atoms. It does not provide a neighbourlist for the ghost atoms. This corresponds to consider_ghost_neighbours=false. And it is fine if only a pair list is needed for whatever comes afterwards. But if triplets are constructed subsequently, the neighbours of ghost atoms are needed. Some atomic representations have three-body terms and they can not be constructed if the ghost atoms do not have neighbours.
AdaptorHalfList
--------------
AdaptorHalfList provides to functionality to reduce a full neighbourlist to a half list. After construction the iteration over all possible pairs in an atomic structure yields each pair only once. The choice for including or not including a pair is based on the atomic indices of the participating atoms. Only pairs where the first (i-atom) index is smaller than the second (j-atom) index are included. Atomic indices are rascal-internal and the user does not have control over atomic indices. This ensures that the condition is a safe check.
AdaptorFullList
-----------------
AdaptorFullList provides the functionality to blow up a half neighbourlist to a full neighbourlist. After construction the iteration over pairs yields each atom at least once as an i-atom. By construction the logic of layering in rascal for pairs is broken. Hence layer information for all pairs (also the ones that previously existed) is reset.

This means that a previously calculated Property is not accessible any more with :class:`~ClusterRef`s of this type. Or it is partly accessible, but meaningless since the layering is destroyed.
AdaptorStrict
--------------
AdaptorStrict provides the functionality to only iterate over neighbours that are actually within a cutoff. Can must be stacked on top of an existing neighbourlist. During calculation of the strict neighbourlist the actual distance between neighbouring atoms has to be calculated. The distance (often used for calculating representations) and direction (used in derivatives) is stored for later possible reuse.

Using this adaptor ensures that later calculation of (possibly very expensive) representations is limited to the atoms which are in each other’s cutoff.
AdaptorMaxOrder
------------------
AdaptorMaxOrder increases and existing neighbourlist by one Order. I.e. a pair list becomes a triplet list, a triplet list a quadruplet list, etc. It must be stacked on an existing neighbourlist.

The functionality is provided by taking a Cluster of a given order and adding all neighbours (pairs) of the constituting atoms as neighbours of this cluster.

(Currently there is a peculiarity. Adding all neighbours results in an ambiguity. The triplet list is not strict any more. This is expected to be sorted out during the implementation of Representations)
AdaptorFilter

AdaptorFilter is a pure virtual which can not be instantiated. An daughter class of it provides a filtered (masked) view on an existing StructureManager or Adaptor.

Any daughter class needs to implement the perform_filtering() function. When used in this way, it could in principle substitute the AdaptorStrict.

Currently only the SpeciesManager uses this class. Within this context the AdaptorFilter is used as a bucket to store a species-sorted view.

Species Manager
===============

Principle
-------------
A SpeciesManager is an Object which provides access to a structure sorted by species. Upon construction the whole tree of combinations is unrolled up to the MaxOrder (i.e. pairs or triplets typically). The use case is e.g. a function (symmetry function, etc.) which is defined for a given atom species combination with types 1-2-1. Access is provided via square bracket operater taking an array of integers, i.e. std:array<int, 3> species_indices{1, 2, 1}. The return type is an iterable filtered view on the structure sorted by species_indices. Internally this is converted to a KeyStandardisation.
Storage concept
---------------
Internally the SpeciesManager stores the splitted species in a container that is a map of KeyStandardisation and AdaptorFilter (please refer to adaptor filter for details). `AdaptorFilter`s serve as a bucket to store atoms/pairs/triplets of a given combinations.

The choice of implementation is based on two ideas #. Speed #. Reuseability

Sorting a given StructureManager with SpeciesManager into specific combinations yields in one iteration over all the possible Clusters, i.e. maximum speed for sorting. If a new combination of species if found during iteration, a new bucket is added to the map and the respective cluster is put into that bucket. If a respective bucket already exists, the cluster is placed in that.

Putting them into containers of AdaptorFilter makes it easy to reuse this class type for other kinds of sorting.

A downside of this choice is that the concept of update_self() is broken. In the current case the AdaptorFilter just serves as a bucket to put a specific type of cluster in. It does not know anything about the rationale why it was given this combination. So the brains here lie with the StructureManager. Hence, then recursive update (if e.g. the underlying structure changes) can only be invoked by the StructureManager who does the actual sorting. Further, if a new structure is provided, it could be that the AdaptorFilter, who was the container for a specific combination of species, does not exist anymore because the neighbourhood changed.
Access
------------
One example for how to access a filtered view on a structure through the SpeciesManager is given in the test. By using a std::array<int, Order> filled with the desired species combination with the operator[], a StructureManager is returned, which only holds this specific combination.

