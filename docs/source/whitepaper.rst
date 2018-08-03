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
It takes care of 

* Holding information about an atomic structure
* Building clusters based on existing structural information
* Provide iterators that can be used to loop over clusters of a given type,
and to access properties associated with clusters of each body order.


