.. _cluster_indices_container:

Cluster indices container
-------------------------

.. Here cluster indices container object explained ignoring number of layers,

.. Reason for cluster indices container

To iterate fast over clusters, a list of cluster indices of a required order is stored in a contiguous memory structure.
The lists of cluster indices for all required orders are stored in the `cluster_indices_container` member variable of the :class:`StructureManager` class.

.. here I am not sure about the exact reasons of performance,
  but I think the performance reasons should be explained in a different section, where also the meta programming and its performance reasons are explained, and then referenced here

Due to reasons of performance, we use :class:`~Property` objects to access and store a list of cluster indices.
For each order one :class:`~Property` object is constructed and stored in the `cluster_indices_container`.
Therefore `cluster_indices_container` is a tuple of :class:`~Property` objects with size of the maximal needed order.
For now we only consider a one layer structure manager, to explain the access on the cluster indices using `cluster_indices_container`.

.. literalinclude:: ../../../src/structure_managers/adaptor_neighbour_list.hh
    :language: c++
    :start-after: adaptor-neighbour-list-get-atom-pair-properties-start
    :end-before: adaptor-neighbour-list-get-atom-pair-properties-end
    :caption: src/structure_managers/adaptor_neighbour_list.hh :class:`~AdaptorNeighbourList`
    :dedent: 4

.. maybe the above quoted code will be made to a distinct file, where this principle is more explained.  then this part can be copied to the new files 
  .. code-block:: c++
  
      // the property which stores an atom index list in its values successively
      // {atom_0_index, atom_1_index, ...}
      auto & atom_index_property{std::get<0>(manager->cluster_indices_container)} 
      // 
      // the property which stores a pair index list in its values as 
      // {(neigbour_0_of_atom_0, atom_0)_index, (neigbour_1_of_atom_0, atom_0)_index, ...,
      //     (neigbour_0_of_atom_1, atom_1)_index, (neigbour_1_of_atom_1, atom_1), ...}
      auto & pair_index_property{std::get<1>(manager->cluster_indices_container)} 

In this way if one wants to iterate through all pairs of atoms for some atom, it can be iterated contiguously through the pair index list giving an the starting index. This can be further generalized to get the cluster indices lists of any order `Order` as seen below. 

.. literalinclude:: ../../../src/structure_managers/adaptor_increase_maxorder.hh
    :language: c++
    :start-after: adaptor-increase-maxorder-get-maxorder-property-start
    :end-before: adaptor-increase-maxorder-get-maxorder-property-end
    :caption: src/structure_managers/adaptor_increase_maxorder.hh :class:`~AdaptorMaxOrder`
    :dedent: 4

.. maybe the above quoted code will be made to a distinct file, where this principle is more explained.  then this part can be copied to the new files 
  .. code-block:: c++
      // the property which stores a cluster index list of order Order in its values as 
      // {(neigbour_0_of_cluster<Order-1>_0, cluster<Order-1>_0)_index,
      //     (neigbour_1_of_cluster<Order-1>_0, cluster<Order-1>_0)_index, ...,
      //     (neigbour_0_of_cluster<Order-1>_1, cluster<Order-1>_1)_index,
      //     (neigbour_1_of_cluster<Order-1>_1, cluster<Order-1>_1), ...}
      // be awared that the cluster_indices_container tuple is numbered starting at 0,
      // therefore to get the property of order Order, we have to access position Order-1
      auto & cluster_index_property{std::get<Order-1>(manager->cluster_indices_container)} 

The definition of the neighbour of a cluster can be anything from the nearest atom in terms of the smallest sum of distance to all atoms in the cluster, up to all atoms in cell plus ghost atoms.

When manager implementations are stacked, the `cluster_indices_container` tuple is copied to the new structure manager. An access request on the cluster indices can be therefore resolved inside the current level structure manager and has not to be forwarded to lower level structure managers. 

Iteration through clusters
**************************

The class :class:`~StructureManager` and the nested classes :class:`~StructureManager::Iterator` and :class:`~StructureManager::ClusterRef` in `structure_manager.hh` are essential to understand the iteration through clusters. We use their nested name as shortcut. The iteration through clusters of order 1 (atoms) is done by iterating over an structure manager as seen below. 

.. code-block:: c++

    for (auto atom : manager) {
      // do something with atom
    }

By using a structure manager as iterator the `begin()` member function of the manager is invoked. This function is inherited from the structure manager definition to each structure manager implementation. In the function a :class:`~Iterator` object is constructed. 

.. literalinclude:: ../../../src/structure_managers/structure_manager.hh
    :language: c++
    :start-after: structure-manager-iterator-start
    :end-before: structure-manager-iterator-end
    :caption: src/structure_managers/structure_manager.hh :class:`~StructureManager`
    :dedent: 4

This iterator object is always constructed of type order 1.

.. literalinclude:: ../../../src/structure_managers/structure_manager.hh
    :language: c++
    :start-after: structure-manager-iterator-types-start
    :end-before: structure-manager-iterator-types-end
    :caption: src/structure_managers/structure_manager.hh :class:`~StructureManager`
    :dedent: 4

The iterator only increases its `index` member variable in each iteration step and returns a pointer of its object. Thereby invoking its `operator*[]` member function, which returns a :class:`~ClusterRef` object (in the example above, this is the `atom` object). This object holds the information about the cluster indices for each layer and the atom indices in its cluster corresponding to the `index` of the iterator.

.. literalinclude:: ../../../src/structure_managers/structure_manager.hh
    :language: c++
    :start-after: iterator-operator-pointer-access-start
    :end-before: iterator-operator-pointer-access-end
    :caption: src/structure_managers/structure_manager.hh in :class:`~Iterator` 
    :dedent: 4

For the construction of this object the iterator gets the property in the `cluster_indices_container` corresponding to the order of the current iteration. `Order-1` is used because the tuple is numbered starting from 0. The property is accessed with the cluster index, which is the sum of an offset and the iterators index. In the case of clusters of order 1 (atoms) the iterator is always constructed with offset 0.

.. literalinclude:: ../../../src/structure_managers/structure_manager.hh
    :language: c++
    :start-after: iterator-get-cluster-index-start
    :end-before: iterator-get-cluster-index-end
    :caption: src/structure_managers/structure_manager.hh in :class:`~Iterator` 
    :dedent: 4

The rest is just a look up in the property values with the cluster index.
When it is looped over a :class:`~ClusterRef` object as with `atom` in

.. code-block:: c++

    for (auto atom : manager) {
      for (auto pair : manager) {
        // do s
      }
    }

a new iterator is created, which then returns new :class:`ClusterRef` objects for each iteration step.
The :class:`ClusterRef` and Iterator objects always store a reference of its creator object.
There is a chain of linked iterators and cluster ref objects when doing nested loops.
This is important to calculate the offsets for clusters.
The offsets of the new Iterator is assigned in the `begin()` function of :class:`ClusterRef` with the indices of all currently nested loops as argument called counters. 

.. literalinclude:: ../../../src/structure_managers/structure_manager.hh
    :language: c++
    :start-after: clusterref-iterator-start
    :end-before: clusterref-iterator-end
    :caption: src/structure_managers/structure_manager.hh in :class:`ClusterRef`` 
    :dedent: 4

The calculation of the offsets is part of the manager implementation. To determine the end of the loop the `size` function is invoked which is also forwards this access request to the manager implementation. With each new nested loop the procedure repeats itself.


Layering
********

Here we introduce how the layers change the structure of the cluster properties.

.. Reason for layering

To only iterate over clusters with certain properties, while still iterating contiguously, an adaptor filter can be used.
Within this kind of adaptor a new cluster index list is created from the lower level list containing only the clusters with the needed property.
To still be able to reference clusters with their original index, the original index together with the new cluster index of the filtered clusters is stored inside the adaptor's cluster indices list property of the corresponding order.
In addition the adaptor filter can store physical properties inside a property which can be accessed with its new cluster index.
To be able to reference properties from lower layers, the cluster indices of the filtered clusters from all lower layer managers are copied up to each current layer manager.
Effectively, each stacked adaptor even copies, copies a subset and expands, or resets the cluster indices list for all orders. For each order a different choice can be made.
If the adaptor filters clusters, it only copies the indices from the filtered clusters and expands them by a new index.
By this procedure the adaptor of the highest level has a copy of all cluster indices with calculated properties. 
Therefore the number of cluster indices per cluster is the same as the manager's number of layers.
This is one the reason for the concept of layers. 

Below it is visualized how cluster indices lists for managers with different numbers of layers is stored.
In this view the order of the list is arbitrary but fixed. The most left cluster index of a manager with a corresponding layer is the new generated cluster index list. In the table CI_i means cluster index of layer i. For a layer a row is stored in the `cluster_indices` member variable of the :class:`~ClusterRef` object when iterating over a structure manager with this amount of layers.

.. TODO maybe do this in code to be more understandable

+-----------+-------------+--------------------+
| One layer | Two layers  |   Three layers     |
+-----------+------+------+------+------+------+
| CI_0      | CI_1 | CI_0 | CI_2 | CI_1 | CI_0 |
+===========+======+======+======+======+======+
|   0       | 0    | 1    | 0    | 1    | 2    |
+-----------+------+------+------+------+------+
| 1         | 1    | 2    |      |      |      |
+-----------+------+------+------+------+------+
| 2         |      |      |      |      |      |
+-----------+------+------+------+------+------+

For example a :class:`~StructureManagerCenters` object with 3 atoms would have for order 0 the following list.

+-----------+
| One layer |
+-----------+
| CI_0      |
+===========+
|   0       |
+-----------+
| 1         |
+-----------+
| 2         |
+-----------+

When iterating over the structure manager, the `cluster_indices` of the :class:`~ClusterRef` objects would a have a size of 1. The first object would contain the atom with index 0 and so on.
Then stacking one filter adaptor `AdaptorSomeFilter1<StructureManagerCenters>` filtering the atom with index 2 results in the cluster_indices_list:

+-------------+
| Two layers  |
+------+------+
| CI_1 | CI_0 |
+======+======+
| 0    | 1    |
+------+------+
| 1    | 2    |
+------+------+

The cluster indices list property stores now a list of two cluster indices and is therefore of layer 2 and uses 2 entries (the property's `NbRow` template argument is 2) to store the cluster indices. 
Stacking another filter adaptor `AdaptorSomeFilter2<AdaptorSomeFilter1<StructureManagerCenters>>` filtering the atom with index 1 results in: 

+--------------------+
|   Three layers     |
+------+------+------+
| CI_2 | CI_1 | CI_0 |
+======+======+======+
| 0    | 1    | 2    |
+------+------+------+

When iterating over this structure manager, only one :class:`~ClusterRef` object would be returned with `cluster_indices` {0,1,2}.

..
  Why no forwarding? (in more detail)
  forward to lower layer, then ClusterRef are returned from lower layer which cannot be used for iterating in higher levels
  Using forwarding with the CRTP would lead to less code flexibility, because it is allowed to stack adaptors in a variety of ways.
