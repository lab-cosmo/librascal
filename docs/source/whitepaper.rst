.. _whitepaper:

Implementation 
=================

Rascal is built as an efficient C++ core, that can be easily accessed 
from Python. Coding philosophy and style conventions are discussed in the
Developer's guide dev_

Here we discuss the general structure of the code, the main objects that
are used, and the rationale for various design choices.

Structures and Clusters
-----------------------

An atomic structure is a collection of atoms, with properties associated 
to it. In the most abstract way possible, properties can be attached to the
entire structure, to atoms, or to groups of atoms (pairs, triplets, ...)
which we refer to as "clusters".
For instance, the energy is a property associated to an entire molecule or
crystal structure, the nuclear chemical shift is a property associated to 
atoms, a distance is a property of a pair of atoms, and so on. 
We refer to the number of atoms involved in a 

Structure Managers
------------------


``StructureManager`` is the main class forming the backbone of Rascal. 
It takes care of:

* Holding information about an atomic structure
* Building clusters based on existing structural information
* Provide iterators that can be used to loop over clusters of a given type,
and to access properties associated with clusters of each body order.

.. image:: ../figures/implementation_structure_manager.svg

In practice the ``StructureManager`` allocates the memory to store
property data, which can be accessed as a heterogeneous map containing
different types of properties. It also provides an interface to access 
clusters in terms of iterators: the internal implementation of the 
neighbor list is left to the implementation, and it only should allow
to retrieve ``ClusterRef`` objects, that act as a light pointer to the 
sequential list of clusters, and contains the indices of the atoms involved
in the cluster. 

``ClusterRef`` objects can be in turn iterated on, giving access to a 
concise notation to iterate over atoms, pairs, triplets, ...

This is best explained with a code snippet

.. code-block:: c++
    
    StructureManager SM(structure_data); // assume this can compute pairs
    SM.update();   // refreshes the list of pairs
    
    energy = SM.get_real(0, "energy") // structure-global property access
    //access to atomic positions, explicit typing
    pos = SM.get_3vec(1, "positions"); 
    // pos = SM.get("positions");  // implicit (slower)
    dist = SM.get_real(2, "distances");
    
    for (auto atom: SM) {
      // iterates over order-1 clusters (atoms)
      ri = pos[atom] //accesses the position of an iterated atom
      
      for (auto pair: atom) {
        // iterates over order-2 clusters (pairs)
        dij = dist[pair]  // accesses the distance between the two atoms
        
        rj = pos[pair]  // the cluster can be used to access properties of 
                        // the last member of the cluster
      }
    }
    

Stacking of StructureManagers
-----------------------------

Structure managers can be build in a modular fashion by stacking them on
top of each other. That is, one can add functionalities to an existing 
manager by creating an adapter, which is just a structure manager that is 
built based on a lower-level manager. 
An example would be selecting only atoms of a given kind, or building 
clusters of high body order based on lower-order clusters. 
This leads to the problem that clusters can be defined at different levels:
for instance, one could have pairs defined up to a cutoff of 5 Angstrom, 
and then create a `AdaptorStrict` manager to select only the pairs within 
3 Angstrom, by stacking it on top of the former. `ClusterRef`s to pairs 
iterated within the adapter should point at this reduced list, but single-atom
clusters do not need to be duplicated, as they are still present and valid
in the parent manager.





