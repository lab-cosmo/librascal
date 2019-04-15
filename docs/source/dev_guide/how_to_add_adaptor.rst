.. _how_to_add_adaptor:


Iteration process
-----------------
The SM/CRef returns a SM::iterator when used for iteration with itself as argument. Then the iterator only increases its index member variable in each iteration ++i by this with operator*[] it returns CR with the access cluster_indices by accessing the property of the underlying SM with the iterators order and current index.



Container_t & container;
the iterator has a container showing the SM/CR constructed him 
How to add a new adaptor
------------------------

Because of the CRTP scheme, there are several member functions which are enforced to be implementey by an adaptor.
Besides these enforced member functions, the implementation depends mostly on the functionality of the adaptor.
Common functionalities like adding calculating a property of clusters are also explained here.

TODO explanation
- traits
- ghost atoms

structure traits which have to be implemented
* Dim: The dimensionality of the clusters in the manager. A value of 3 means that 3 dimensional data for the clusters are supported.
* LayerByOrder: It is a sequence of integers (unsigned integers of type size_t) where the position in the sequence is represents the order and the integer at this position expresses the number of layers with this orders. The number of layers for an order is counted according to zero-based numbering. For example the sequence `std::index_sequence<0, 5, 2>` means that there is 1 layer for order 1, 6 layers for order 2 and 3 layers for order 3. There are helper functions in `CRK` `LayerExteended` for extending the size of the sequence when increasing the maximal order and `LayerIncreaser` for incrementing for each order the number of layers. 

used traits
    constexpr static int Dim{3}; // always required
    constexpr static size_t MaxOrder{2};
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{false};
    using LayerByOrder = std::index_sequence<0, 0>; // always required


This is for structure managers, adaptors also have to implement these ones, but can use the underlying function of the manager
Due to CRTP a derived class offers implementations of the interface functions but also new functionalities.

******
If one of the enforced member functions does not fit to the functionality of the structure manager like in the case of the `get_cluster_neighbour` function for an structure manager with order 1, there has still to be an implementation catching this case as for example in structure_manager_centers.cc
    template <size_t Order, size_t Layer>
    inline int
    get_cluster_neighbour(const ClusterRefKey<Order, Layer> & /*cluster*/,
                          size_t) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms.");
      return 0;
    }

functions required implementation:
get_size in interface size()
get_size_with_ghosts
get_nb_clusters
get_cluster_size
get_cluster_neighbour // only for apator of order 2 important // important for iteration process // implemented by the sm
get_position
get_atom_type
get_offset_impl

unlike virtual function not to implementat these functions will not result in an compiler error, but a user should always be able to access this function trough the sm interface and if an interface function does not make sense for a certain adaptor/sm it should throw a warning and a default value as in smc with get_cluster_size function

nonenforced but used functions
update

.. _add-property-to-adaptor:

How to add a property to an adaptor
****************************

The general usage of a property in an adaptor is explained here on the example of the distance property in ``AdaptorStrict``. As a reminder the ``AdaptorStrict`` calculates for pairs of atoms within a cutoff distance the distance. Therefore it can only be stacked on top of a structure manager of at least order 2. Since an adaptor is always assigned to a layer with its construction, a helper type definition ``Property_t`` is used to automatically assign the layer number to the property in class ``StructureManager``.

.. literalinclude:: ../../../src/structure_managers/structure_manager.hh
    :language: c++
    :start-after: property_t-typedef-start
    :end-before: property_t-typedef-end
    :caption: src/structure_managers/structure_manager.hh
    :dedent: 4
   
Then in the adaptor we specify the property value type `T=double`, because we want to express real numbers with it, and the non-type template parameter `Order=2`, because pairs of atoms are required to calculate distances, and for the other parameters we use the default values as seen above, because a scalar is 1x1 dimensional value. Further, we offer an getter function to access the property from a structure manager given a ``ClusterRef`` object. The code describing the beforementioned can be found in the class ``AdaptorStrict`` as seen below.

.. literalinclude:: ../../../src/structure_managers/adaptor_strict.hh
    :language: c++
    :start-after: property-distance-getter-start
    :end-before: property-distance-getter-end
    :caption: src/structure_managers/adaptor_strict.hh
    :dedent: 4

The distance property is declared in the ``AdaptorStrict`` class as seen below.

.. literalinclude:: ../../../src/structure_managers/adaptor_strict.hh
    :language: c++
    :start-after: property-distance-declaration-start
    :end-before: property-distance-declaration-end
    :caption: src/structure_managers/adaptor_strict.hh
    :dedent: 4

It is constructed in the constructor of ``AdaptorStrict`` with a pointer to the underlying structure manager.

.. literalinclude:: ../../../src/structure_managers/adaptor_strict.hh
    :language: c++
    :start-after: property-distance-construction-start
    :end-before: property-distance-construction-end
    :caption: src/structure_managers/adaptor_strict.hh
    :dedent: 2

In the update function of ``AdaptorStrict`` the actual property values are calculated and assigned to the property. In simplified form this can be written as 

.. code-block:: c++
    :caption: pseudocode of the calculation of the distance property 
    
    template <class ManagerImplementation>
    void AdaptorStrict<ManagerImplementation>::update() {
      for auto atom : this->structure_manager { // simplified
        for auto pair : atom {
          // the distance value
          double distance{calculate_distance(pair)}; // simplified 
          if (distance_value <= cutoff) { // simplified
            // distance value is pushed into the distance property object
            this->distance.push_back(distance);
          }
        }
      }
    }
