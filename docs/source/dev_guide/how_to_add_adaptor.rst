.. _how_to_add_adaptor:

.. BAUSTELLE
  Updating
  ---
  Intialisation, upper level manager forwards update request to lower level up to underlying structure manager, then the structure manager updates the stack of adaptors from lower level to higher level.
  
  
  As in the case of the AdaptorStrict only clusters of order 2 are filtered, but clusters of order 1 or higher order than 2 are ignored.
  In these cases 
  To allow the access of lower level adaptors calculated properties, it is needed to forward enough information from the lower to the upper layers, when inherited.
  To access a property from the lower levels, the property itself has to be known, as well as the cluster indices at that layer.
  For example in a scenario where we have 5 layers for some order and we want to access a property which has been calculated in layer 1, there 
  AF<AF<AF<>> manager
  a property has been calculated in layer 1 for certain scope of clusters, but in layer 5 more clusters have been filtered out with like it can be done with ``AdaptorFilter``. 
  
  , higher level adaptors need to be able to access data from lower level adaptors.
  
  in every way it will make sense from the functionality
  classes with new functionalities can be stacked
  For each order layer pair there is a list of cluster indices.
  
  Why need number of layers: A property could be calculated in one of the lower layers, to access this property, we require the indices of this layer.
  TODOlater Currently this does not work, because properties in lower layers cannot be accessed, but it should work in future.
      -
  is the 2 from layer 2 just copied    
  
  
  It is differed between adaptors which extend the functionality of an order (like adding ghost atoms, initializing a new order), which requires a new sorting of underlying  cluster indices lists.
  
  and adaptors optionally reduning the access scope on the clusters, filter through them,
  which can be implemented as a new index list indexing the filtered atoms with their index in the underlying list.
  ..example look at AdaptorStrict
  
      auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
      auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};
  
  Why do we copy indices even if no filtering happens? Because of the CRTP inheritance the upper layer needs to know the indices of the lower layers. Even it is copied from the lower layer or the request needs to be forwarded to the lower layers. Usually forwarding results in more complex behaviour from the regarding the implementation and the computational ressources, since in each forwarding step it has to be checked if this adaptor is responsible to resolve the request.
  
  
  which filter from an order list, filtering through clusters with certain criteria and calculate properties of these filtered clusters while filtering is optional.
  type A is resets the layer count always 0 at the specific orders
  type B adds a layer to all existing layers
  
  * LayerByOrder: It is a sequence of integers (unsigned integers of type size_t) where the position in the sequence is represents the order and the integer at this position expresses the number of layers with this orders. The number of layers for an order is counted according to zero-based numbering. For example the sequence `std::index_sequence<0, 5, 2>` means that there is 1 layer for order 1, 6 layers for order 2 and 3 layers for order 3. There are helper functions in `CRK` `LayerExtender` for extending the size of the sequence when increasing the maximal order and `LayerIncreaser` for incrementing for each order the number of layers. 

How to add a new adaptor
------------------------

.. Because of the CRTP scheme, there are several member functions which are enforced to be implementation by an adaptor.
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

  It is important to differ between neigbour0_of_atom0 as pair index as it is stored in the above property, and as atom index as it is stored in the :class:`~AdaptorNeigbourList` in its `neighbours` member variable. The pair indices can be completely independent from the atom indices, while in the other case the indices have to correspond to the indices in the property of order 1. For that reason ClusterRef has the the atom indices as additional member variable.

.. _add-property-to-adaptor:

.. TODO this section should be somewhere else

Adaptor filter and extender
***************************

When a new layer is added on top of some order, then the cluster indices are not only copied, but they are expanded by an cluster index, which indexes the new layer.
Therefore adaptors which filter clusters always increase the number of layers for each existing order. These kind of adaptors are called adaptor filters.

On the other hand adaptors which add new clusters reset the layer only on the orders where they are applied on.
The reason for this is that the cluster indices in the ground layer have to be changed, and the filters in upper layers are not recalculated to prevent a forwarding of such a request. As an example the AdaptorNeigbourList calculates atoms and their neigbours and saves their indices in the cluster indices of order 2. This action resets the order 2 layers.
In addition ghost atoms are needed in the neighbour calculation to calculate all pairs, therefore the ghost atoms also become part of the atom indices and have to be reseted to layer 0.
These kind of adaptors are called adaptor extender in the sense that the adaptor extends the functionality of certain clusters. For example the AdatporNeigbourList extends the clusters of order 1 by ghost atoms. 




How to add a property to an adaptor
***********************************

This section has been commented out, because it was not correct and has to be reviewed

.. commented out
    The general usage of a property in an adaptor is explained here on the example of the distance property in ``AdaptorStrict``. As a reminder the ``AdaptorStrict`` calculates for pairs of atoms within a cutoff distance the distance. Therefore it can only be stacked on top of a structure manager of at least order 2. Since an adaptor is always assigned to a layer with its construction, a helper type definition ``Property_t`` is used to automatically assign the layer number to the property in class ``StructureManager``.
    
    .. literalinclude:: ../../../src/structure_managers/structure_manager.hh
        :language: c++
        :start-after: property_t-typedef-start
        :end-before: property_t-typedef-end
        :caption: src/structure_managers/structure_manager.hh
        :dedent: 4
       
    Then in the adaptor we specify the property value type `T=double`, because we want to express real numbers with it, and the non-type template parameter `Order=2`, because pairs of atoms are required to calculate distances, and for the other parameters we use the default values as seen above, because a scalar is 1x1 dimensional value. Further, we offer an getter function to access the property from a structure manager given a ``ClusterRef`` object. The code corresponding to the beforementioned description can be found in the class ``AdaptorStrict`` as seen below.
    
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
