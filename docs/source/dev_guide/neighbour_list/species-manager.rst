.. _species-manager:

Species Manager
~~~~~~~~~~~~~~~

.. contents::
   :local:


Principle
*********
A :class:`~SpeciesManager` is an Object which provides access to a structure sorted by species. Upon construction the whole tree of combinations is unrolled up to the ``MaxOrder`` (i.e. pairs or triplets typically). The use case is e.g. a function (symmetry function, etc.) which is defined for a given atom species combination with types 1-2-1.
Access is provided via square bracket operater taking an array of integers, i.e. ``std:array<int, 3> species_indices{1, 2, 1}``. The return type is an iterable filtered view on the structure sorted by ``species_indices``. Internally this is converted to a :class:`~KeyStandardisation`.

Storage concept
***************
Internally the :class:`~SpeciesManager` stores the splitted species in a container that is a map of :class:`~KeyStandardisation` and :class:`~AdaptorFilter` (please refer to :ref:`adaptor filter` for details).
`AdaptorFilter`s serve as a bucket to store atoms/pairs/triplets of a given combinations.

The choice of implementation is based on two ideas
#. Speed
#. Reuseability

Sorting a given :class:`~StructureManager` with :class:`~SpeciesManager` into specific combinations yields in one iteration over all the possible Clusters, i.e. maximum speed for sorting. If a new combination of species if found during iteration, a new `bucket` is added to the map and the respective cluster is put into that `bucket`. If a respective `bucket` already exists, the cluster is placed in that.

Putting them into containers of :class:`~AdaptorFilter` makes it easy to reuse this class type for other kinds of sorting.

A downside of this choice is that the concept of :meth:`~update_self()` is broken. In the current case the :class:`~AdaptorFilter` just serves as a bucket to put a specific type of cluster in. It does not know anything about the rationale why it was given this combination. So the `brains` here lie with the :class:`~StructureManager`. Hence, then recursive update (if e.g. the underlying structure changes) can only be invoked by the :class:`~StructureManager` who does the actual sorting.
Further, if a new structure is provided, it could be that the :class:`~AdaptorFilter`, who was the container for a specific combination of species, does not exist anymore because the neighbourhood changed.

Access
******
One example for how to access a filtered view on a structure through the :class:`~SpeciesManager` is given in the test. By using a ``std::array<int, Order>`` filled with the desired species combination with the ``operator[]``, a :class:`~StructureManager` is returned, which only holds this specific combination.
