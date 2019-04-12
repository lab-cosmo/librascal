.. _species-manager:

Species Manager
~~~~~~~~~~~~~~~

.. contents::
   :local:


Principle
*********
A `SpeciesManager` is an Object which provides access to a structure sorted by species. Upon construction the whole tree of combinations is unrolled up to the `MaxOrder` (i.e. pairs or triplets typically). The use case is e.g. a function (symmetry function, etc.) which is defined for a given atom species combination with types 1-2-1.
Access is provided via square bracket operater tacking an array of integers, i.e. `std:array<int, 3> species_indices{1, 2, 1}`. The return type is an iterable filtered view on the structure sorted by `species_indices`.

Storage concept
***************
Internally the `SpeciesManager` stores the splitted species in a container that is a map of `KeyStandardisation` and `AdaptorFilter` (please refer to :ref:`adaptor filter` for details).
`AdaptorFilter`s serve as a bucket to store atoms/pairs/triplets of a given combinations.

The choice
