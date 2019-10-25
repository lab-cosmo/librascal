/**
 * @file   adaptor_filter.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   23 Oct 2018
 *
 * @brief An adaptor that provides a filtered (masked) view
 *        on an existing structure manager.
 *
 * Copyright  2018 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_FILTER_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_FILTER_HH_

#include "rascal_utility.hh"
#include "structure_managers/cluster_ref_key.hh"
#include "structure_managers/filter_base.hh"
#include "structure_managers/structure_manager.hh"

#include <type_traits>

namespace rascal {

  /**
   * Forward declaration for traits
   */
  template <class ManagerImplementation, size_t MaxOrder>
  class AdaptorFilter;

  /**
   * Specialisation of traits for increase <code>MaxOrder</code> adaptor
   */
  template <class ManagerImplementation, size_t MaxOrder_>
  struct StructureManager_traits<
      AdaptorFilter<ManagerImplementation, MaxOrder_>> {
    using parent_traits = StructureManager_traits<ManagerImplementation>;
    constexpr static AdaptorTraits::Strict Strict{parent_traits::Strict};
    constexpr static bool HasDistances{parent_traits::HasDistances};
    constexpr static bool HasDirectionVectors{
        parent_traits::HasDirectionVectors};
    constexpr static bool HasCenterPair{parent_traits::HasCenterPair};
    constexpr static int Dim{parent_traits::Dim};
    constexpr static int StackLevel{parent_traits::StackLevel + 1};
    //! New MaxOrder upon construction!
    constexpr static size_t MaxOrder{MaxOrder_};
    //! New Layer
    //! TODO: Is this the correct way to initialize the increased order?
    using LayerByOrder =
        typename LayerIncreaser<MaxOrder_,
                                typename parent_traits::LayerByOrder>::type;
  };

  namespace internal {

    /**
     * When a cluster is added, it is not generally known whether the parent
     * cluster (e.g, a pair in the case of a triplet, or the atom in the case of
     * the pair) has already been added. This structure checks this recursively.
     */
    template <size_t Order, class ManagerImplementation, size_t MaxOrder>
    struct ClusterAdder {
      using Manager_t = AdaptorFilter<ManagerImplementation, MaxOrder>;
      using ClusterRef_t =
          typename Manager_t::template InputClusterRef_t<Order>;
      static void add_parent(Manager_t & manager,
                             const ClusterRef_t & cluster) {
        // e.g., the pair (i,j) for a triplet (i,j,k), or the atom i
        // for a pair (i,j)
        const auto & parent_cluster{cluster.get_iterator().get_container()};
        if (not manager.has_cluster(parent_cluster)) {
          manager.add_cluster(parent_cluster);
        }
      }
    };

    /**
     * Recursion end, where nothing has to be done
     */
    template <class ManagerImplementation, size_t MaxOrder>
    struct ClusterAdder<1, ManagerImplementation, MaxOrder> {
      static constexpr size_t Order{1};
      using Manager_t = AdaptorFilter<ManagerImplementation, MaxOrder>;
      using ClusterRef_t =
          typename Manager_t::template InputClusterRef_t<Order>;
      static void add_parent(Manager_t & /*ignored manager*/,
                             const ClusterRef_t & /*ignored atom*/) {}
    };

  }  // namespace internal

  /**
   * Simple Adaptor that can be stacked on top of any StructureManager of lower
   * or equal MaxOrder, and present a filtered view on it. The filtering can be
   * performed based on arbitrary criteria using the add_cluster() method.
   */
  template <class ManagerImplementation, size_t MaxOrder>
  class AdaptorFilter
      : public FilterBase,
        public StructureManager<AdaptorFilter<ManagerImplementation, MaxOrder>>,
        public std::enable_shared_from_this<
            AdaptorFilter<ManagerImplementation, MaxOrder>> {
   public:
    template <size_t Order_, class ManagerImplementation_, size_t MaxOrder_>
    friend struct internal::ClusterAdder;

    using Manager_t = AdaptorFilter<ManagerImplementation, MaxOrder>;

    using Parent = StructureManager<Manager_t>;
    using ParentBase = FilterBase;
    using ManagerImplementation_t = ManagerImplementation;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using traits = StructureManager_traits<AdaptorFilter>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    using Vector_ref = typename Parent::Vector_ref;
    template <size_t Order>
    using InputClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;
    using Hypers_t = typename Parent::Hypers_t;

    static_assert(traits::MaxOrder <= ManagerImplementation::traits::MaxOrder,
                  "can only present view on existing clusters");

    //! Default constructor
    AdaptorFilter() = delete;

    //! constructor underlying manager
    explicit AdaptorFilter(ImplementationPtr_t manager) : manager{manager} {
      this->reset_initial_state();
    }

    //! Copy constructor
    AdaptorFilter(const AdaptorFilter & other) = delete;

    //! Move constructor
    AdaptorFilter(AdaptorFilter && other) = default;

    //! Destructor
    virtual ~AdaptorFilter() = default;

    //! Copy assignment operator
    AdaptorFilter & operator=(const AdaptorFilter & other) = delete;

    //! Move assignment operator
    AdaptorFilter & operator=(AdaptorFilter && other) = default;

    /**
     * clears the state fully without deallocating any memory. Needs to be
     * called *before* adding clusters (i.e., also at the beginning of every
     * update)
     */
    void reset_initial_state();

    //! updates the underlying adaptor
    void update_self() {
      this->reset_initial_state();
      this->perform_filtering();
    }

    /**
     * perform the actual filtering work, or send a signal to whoever performs
     * this work. Needs to be implemented in the daughter class
     */

    virtual void perform_filtering() = 0;

    /**
     * return the number of 'neighbours' (i.e., number of pairs for an atom,
     * number of triplets for a pair, etc) of a given order.
     */
    size_t get_nb_clusters(int order) const {
      return this->atom_tag_list[order - 1].size();
    }

    /**
     * return the number of atoms
     */
    size_t get_size() const { return this->get_nb_clusters(1); }

    /**
     * return the position of a given atom
     */
    Vector_ref get_position(int index) {
      return this->manager->get_position(index);
    }

    //! returns the number of atoms
    size_t get_size_with_ghosts() const {
      return this->manager->get_size_with_ghosts();
    }

    //! returns the distance between atoms in a given pair
    template <size_t Order, size_t Layer,
              bool HasDistances = traits::HasDistances>
    const std::enable_if_t<HasDistances, double>
    get_distance(const ClusterRefKey<Order, Layer> & pair) const {
      static_assert(HasDistances == traits::HasDistances,
                    "HasDistances is used for SFINAE, please don't specify it");
      return this->manager->get_distance(pair);
    }

    /**
     * return direction vector
     */
    template <size_t Order, size_t Layer,
              bool HasDistances = traits::HasDistances>
    std::enable_if_t<HasDistances, Vector_ref>
    get_direction_vector(const ClusterRefKey<Order, Layer> & pair) const {
      static_assert(HasDistances == traits::HasDistances,
                    "HasDistances is used for SFINAE, please don't specify it");
      return this->manager->get_direction_vector(pair);
    }

    //! get atom_tag of index-th neighbour of this cluster
    template <size_t Order, size_t Layer>
    int get_neighbour_atom_tag(const ClusterRefKey<Order, Layer> & cluster,
                               int index) const {
      static_assert(Order <= traits::MaxOrder - 1,
                    "this implementation only handles upto traits::MaxOrder");
      auto && offset = this->offsets[Order][cluster.get_cluster_index(Layer)];
      return this->atom_tag_list[Order][offset + index];
    }

    //! get atom_tag of the index-th atom in manager
    int get_neighbour_atom_tag(const Parent & /*parent*/, size_t index) const {
      return this->atom_tag_list[0][index];
    }

    //! return atom type
    int get_atom_type(const AtomRef_t & atom) const {
      // careful, atom refers to our local index, for the manager, we need its
      // index:
      auto && original_atom{this->atom_tag_list[0][atom.get_index()]};
      return this->manager->get_atom_type(original_atom);
    }

    //! Returns atom type given an atom tag
    int get_atom_type(int atom_id) const {
      auto && type{this->manager->get_atom_type(atom_id)};
      return type;
    }

    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template <size_t Order>
    size_t get_offset_impl(const std::array<size_t, Order> & counters) const {
      return this->offsets[Order][counters.back()];
    }

    //! return the number of neighbours of a given atom
    template <size_t Order, size_t Layer>
    size_t
    get_cluster_size_impl(const ClusterRefKey<Order, Layer> & cluster) const {
      static_assert(Order <= traits::MaxOrder - 1,
                    "Order exceeds maxorder for this filter.");
      constexpr auto nb_neigh_layer{
          compute_cluster_layer<Order>(typename traits::LayerByOrder{})};
      return this->nb_neigh[Order][cluster.get_cluster_index(nb_neigh_layer)];
    }

    /**
     * add a cluster to the filter. This is the main use of this
     * class. You can iterate over the StructureManager you're
     * filtering, and add the atoms/pairs/triplet... you wish to
     * retain in the filtered version using this method
     */
    template <size_t Order>
    void add_cluster(const InputClusterRef_t<Order> & cluster);

   protected:
    /**
     * check whether a cluster has been added already (Only check the
     * last inserted cluster, relying on the underlying iteration
     * order
     */
    template <size_t Order>
    bool has_cluster(const InputClusterRef_t<Order> & cluster);

    /**
     * main function during construction of the filtered view
     * @param cluster last atom of cluster is added to the filter
     * @param Order select whether it is an i-atom (order=1), j-atom (order=2),
     * or ...
     */
    template <size_t Order>
    void add_atom(const InputClusterRef_t<Order> & cluster) {
      const auto & atom_tag{cluster.back()};
      static_assert(Order - 1 <= traits::MaxOrder,
                    "you can only add neighbours to the n-th degree defined by "
                    "MaxOrder of the underlying manager");

      // add new atom at this Order
      this->atom_tag_list[Order - 1].push_back(atom_tag);
      // count that this atom is a new neighbour
      this->nb_neigh[Order - 1].back()++;
      this->offsets[Order - 1].back()++;

      for (auto i{Order}; i < traits::MaxOrder; ++i) {
        // make sure that this atom starts with zero lower-Order neighbours
        this->nb_neigh[i].push_back(0);
        // update the offsets
        this->offsets[i].push_back(this->offsets[i].back() +
                                   this->nb_neigh[i].back());
      }
    }

    //! underlying manager to be filtered
    ImplementationPtr_t manager;

    /**
     * store atom tags per order,i.e.
     *   - atom_tag_list[0] lists all i-atoms
     *   - atom_tag_list[1] lists all j-atoms
     *   - atom_tag_list[2] lists all k-atoms
     *   - etc
     */
    std::array<std::vector<int>, traits::MaxOrder> atom_tag_list{};
    /**
     * store the number of j-atoms for every i-atom (nb_neigh[1]), the number of
     * k-atoms for every j-atom (nb_neigh[2]), etc
     */
    std::array<std::vector<size_t>, traits::MaxOrder> nb_neigh{};
    /**
     * store the offsets from where the nb_neigh can be counted
     */
    std::array<std::vector<size_t>, traits::MaxOrder> offsets{};

   private:
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation, size_t MaxOrder>
  template <size_t Order>
  void AdaptorFilter<ManagerImplementation, MaxOrder>::add_cluster(
      const InputClusterRef_t<Order> & cluster) {
    // The following calls this method recursively on the parent
    // clusters of this cluster, until this method is called with an
    // Order=1 cluster (atom), in which case add_parent does nothing.
    internal::ClusterAdder<Order, ManagerImplementation, MaxOrder>::add_parent(
        *this, cluster);
    this->add_atom(cluster);

    /**
     * Add new Layer for clusters of size Order
     */
    constexpr auto ClusterLayer{
        compute_cluster_layer<Order>(typename traits::LayerByOrder{})};

    Eigen::Matrix<size_t, ClusterLayer + 1, 1> indices{};
    indices.template head<ClusterLayer>() = cluster.get_cluster_indices();
    indices(ClusterLayer) = this->atom_tag_list[Order - 1].size() - 1;
    auto & indices_container{
        std::get<Order - 1>(this->cluster_indices_container)};
    indices_container.push_back(indices);
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation, size_t MaxOrder>
  template <size_t Order>
  bool AdaptorFilter<ManagerImplementation, MaxOrder>::has_cluster(
      const InputClusterRef_t<Order> & cluster) {
    constexpr auto Layer{InputClusterRef_t<Order>::cluster_layer()};

    auto & indices_container{
        std::get<Order - 1>(this->cluster_indices_container)};
    if (indices_container.size() == 0) {
      return false;
    }
    auto && last_cluster_index{indices_container.back()(Layer)};

    return last_cluster_index == cluster.get_cluster_index(Layer);
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation, size_t MaxOrder>
  void AdaptorFilter<ManagerImplementation, MaxOrder>::reset_initial_state() {
    //! Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());
    for (size_t i{0}; i < MaxOrder; ++i) {
      this->atom_tag_list[i].clear();
      this->nb_neigh[i].clear();
      this->offsets[i].clear();
    }

    this->nb_neigh[0].push_back(0);
    for (auto & vector : this->offsets) {
      vector.push_back(0);
    }
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_FILTER_HH_
