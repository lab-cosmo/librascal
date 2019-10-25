/**
 * @file   adaptor_center_contribution.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   July 2019
 *
 * @brief Adaptor that adds the central atom to the iteration over pairs
 *
 * Copyright  2019 Markus Stricker, Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_CONTRIBUTION_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_CONTRIBUTION_HH_

#include "rascal_utility.hh"
#include "structure_managers/property.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/updateable_base.hh"

namespace rascal {
  /*
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorCenterContribution;

  /*
   * specialisation of traits for strict adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<
      AdaptorCenterContribution<ManagerImplementation>> {
    using parent_traits = StructureManager_traits<ManagerImplementation>;
    constexpr static AdaptorTraits::Strict Strict{parent_traits::Strict};
    constexpr static bool HasDistances{parent_traits::HasDistances};
    constexpr static bool HasDirectionVectors{
        parent_traits::HasDirectionVectors};
    constexpr static bool HasCenterPair{true};
    constexpr static int Dim{parent_traits::Dim};
    constexpr static int StackLevel{parent_traits::StackLevel + 1};
    constexpr static size_t MaxOrder{parent_traits::MaxOrder};
    // Reset the Layer for the pair since the underlying iteration pool is
    // increased in size so all the data belongs to the adaptor
    using LayerByOrder = std::index_sequence<
        ManagerImplementation::template cluster_layer_from_order<1>() + 1, 0>;
  };

  /**
   * Adaptor that adds the central atom to the iteration over pairs and the
   * use of the get_atom_ii(). The ii-pair is added at the begining of
   * every block of data associated to neighbors.
   * The default iteration pattern does not include the ii-pair so use
   * with_self_pair() to include it. The properties associated with this
   * manager will have an entry for the ii-pair.
   *
   * It should be stacked after an adaptor providing a linked cell neighbor
   * list, but before any adaptor computing pair related properties as
   * `AdaptoStrict`. This adaptor changes the neighbour list, thus it would
   * make any map of the form pair -> property invalid as it is the case for the
   * distance and direction vector in the `AdaptorStrict` class. The same goes
   * for the `AdaptorMaxOrder` and adaptors computing properties for clusters of
   * Order > 2.
   */
  template <class ManagerImplementation>
  class AdaptorCenterContribution
      : public StructureManager<
            AdaptorCenterContribution<ManagerImplementation>>,
        public std::enable_shared_from_this<
            AdaptorCenterContribution<ManagerImplementation>> {
   public:
    using Manager_t = AdaptorCenterContribution<ManagerImplementation>;
    using Parent = StructureManager<Manager_t>;
    using ManagerImplementation_t = ManagerImplementation;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using traits = StructureManager_traits<Manager_t>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;

    static_assert(traits::MaxOrder == 2,
                  "ManagerImlementation needs to handle pairs");
    static_assert(traits::Strict == AdaptorTraits::Strict::no,
                  "AdaptorCenterContribution does not work on strict NL.");
    constexpr static auto AtomLayer{
        Manager_t::template cluster_layer_from_order<1>()};
    constexpr static auto PairLayer{
        Manager_t::template cluster_layer_from_order<2>()};

    //! Default constructor
    AdaptorCenterContribution() = delete;

    explicit AdaptorCenterContribution(ImplementationPtr_t manager);

    AdaptorCenterContribution(ImplementationPtr_t manager,
                              const Hypers_t & /* adaptor_hypers*/)
        : AdaptorCenterContribution(manager) {}

    //! Copy constructor
    AdaptorCenterContribution(const AdaptorCenterContribution & other) = delete;

    //! Move constructor
    AdaptorCenterContribution(AdaptorCenterContribution && other) = default;

    //! Destructor
    virtual ~AdaptorCenterContribution() = default;

    //! Copy assignment operator
    AdaptorCenterContribution &
    operator=(const AdaptorCenterContribution & other) = delete;

    //! Move assignment operator
    AdaptorCenterContribution &
    operator=(AdaptorCenterContribution && other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    void update_self();

    //! update the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments);

    //! returns the linked cell cutoff for the adaptor
    double get_cutoff() const { return this->manager->get_cutoff(); }

    size_t get_nb_clusters(int order) const {
      return this->atom_tag_list[order - 1].size();
    }

    bool get_consider_ghost_neighbours() const {
      return this->manager->get_consider_ghost_neighbours();
    }

    size_t get_size() const { return this->manager->get_size(); }

    size_t get_size_with_ghosts() const { return this->get_nb_clusters(1); }

    Vector_ref get_position(int index) {
      return this->manager->get_position(index);
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

    /**
     * Since the cluster indices of order 1 are only copied in this filter we
     * can safely use the before-computed list from the previous manager,
     * since they are still valid for access.
     */
    size_t get_atom_index(const int atom_tag) const {
      return this->manager->get_atom_index(atom_tag);
    }

    //! return atom type
    int get_atom_type(const AtomRef_t & atom) const {
      // careful, atom refers to our local index, for the manager, we need its
      // index:
      auto && original_atom{this->atom_tag_list[0][atom.get_index()]};
      return this->manager->get_atom_type(original_atom);
    }

    //! Returns a constant atom type given an atom tag
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
      static_assert(Order < traits::MaxOrder,
                    "Calling this function with the wrong order cluster");
      return this->offsets[Order][counters.back()];
    }

    // TODO(felix) should be removed because not used ?
    template <size_t Order>
    int get_neighbour_atom_tag(const size_t neighbour_index) const {
      static_assert(Order < traits::MaxOrder,
                    "Calling this function with the wrong order cluster");
      return this->manager->get_neighbour_atom_tag(neighbour_index);
    }

    //! return the number of neighbours of a given atom
    template <size_t Order, size_t Layer>
    size_t
    get_cluster_size_impl(const ClusterRefKey<Order, Layer> & cluster) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms and pairs");
      constexpr auto nb_neigh_layer{
          compute_cluster_layer<Order>(typename traits::LayerByOrder{})};
      return this->nb_neigh[Order][cluster.get_cluster_index(nb_neigh_layer)];
    }

    //! Get the manager used to build the instance
    ImplementationPtr_t get_previous_manager() {
      return this->manager->get_shared_ptr();
    }

    // TODO(felix) remove ? or fuse with the function below ?
    /**
     * @return a list of atom tags for clusters of Order 1 aka centers.
     */
    const std::vector<int> get_manager_atom_tag_list() {
      return this->atom_tag_list[0];
    }
    // TODO(felix) rename ?
    /**
     * @return a list of atom tags for clusters of Order 2 aka pairs.
     */
    const std::vector<int> get_neighbours_atom_tag() {
      return this->atom_tag_list[1];
    }

   protected:
    /**
     * main function during construction of a neighbourlist.
     * @param atom_tag the atom to add to the list
     * @tparam Order select whether it is an i-atom (order=1), j-atom (order=2),
     *         etc.
     */
    template <size_t Order>
    void add_atom(int atom_tag) {
      static_assert(Order <= traits::MaxOrder,
                    "you can only add neighbours to the n-th degree defined by "
                    "MaxOrder of the underlying manager");

      // add new atom at this Order
      this->atom_tag_list[Order].push_back(atom_tag);
      // count that this atom is a new neighbour
      this->nb_neigh[Order].back()++;
      this->offsets[Order].back()++;

      for (auto i{Order + 1}; i < traits::MaxOrder; ++i) {
        // make sure that this atom starts with zero lower-Order neighbours
        this->nb_neigh[i].push_back(0);
        // update the offsets
        this->offsets[i].push_back(this->offsets[i].back() +
                                   this->nb_neigh[i].back());
      }
    }

    template <size_t Order, size_t Layer>
    void add_atom(const ClusterRefKey<Order, Layer> & cluster) {
      this->template add_atom<Order - 1>(cluster.back());
    }

    ImplementationPtr_t manager;

    /**
     * For clusters of the form (i,j)
     * store atom tags per order,i.e.
     *   - atom_tag_list[0] lists all i-atoms
     *   - atom_tag_list[1] lists all j-atoms
     */
    std::array<std::vector<int>, traits::MaxOrder> atom_tag_list;
    std::vector<size_t> neighbours_cluster_index;
    /**
     * store the number of j-atoms for every i-atom (nb_neigh[1]), the number of
     * k-atoms for every j-atom (nb_neigh[2]), etc
     */
    std::array<std::vector<size_t>, traits::MaxOrder> nb_neigh;
    /**
     * store the offsets from where the nb_neigh can be counted
     */
    std::array<std::vector<size_t>, traits::MaxOrder> offsets;
  };

  /*--------------------------------------------------------------------------*/
  template <class ManagerImplementation>
  AdaptorCenterContribution<ManagerImplementation>::AdaptorCenterContribution(
      std::shared_ptr<ManagerImplementation> manager)
      : manager{std::move(manager)}, atom_tag_list{},
        neighbours_cluster_index{}, nb_neigh{}, offsets{} {}

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class... Args>
  void AdaptorCenterContribution<ManagerImplementation>::update(
      Args &&... arguments) {
    this->manager->update(std::forward<Args>(arguments)...);
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorCenterContribution<ManagerImplementation>::update_self() {
    //! Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());
    //! initialise the neighbourlist
    for (size_t i{0}; i < traits::MaxOrder; ++i) {
      this->atom_tag_list[i].clear();
      this->nb_neigh[i].clear();
      this->offsets[i].clear();
    }

    this->nb_neigh[0].push_back(0);
    for (auto & vector : this->offsets) {
      vector.push_back(0);
    }

    using AtomIndex_t = typename ClusterRefKey<2, PairLayer>::AtomIndex_t;
    using IndexConstArray =
        typename ClusterRefKey<2, PairLayer>::IndexConstArray;

    // fill the list, at least pairs are mandatory for this to work
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    size_t pair_counter{0};
    // depending on the underlying neighbourlist, the proxy `.with_ghosts()` is
    // either actually with ghosts, or only returns the number of centers.
    for (auto atom : this->manager.get()->with_ghosts()) {
      this->add_atom(atom);
      /*
       * Add new layer for atoms (see LayerByOrder for
       * possible optimisation).
       */
      Eigen::Matrix<size_t, AtomLayer + 1, 1> indices;

      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer - 1);
      atom_cluster_indices.push_back(indices);

      auto && atom_tag = atom.get_atom_tag();
      AtomIndex_t self_atom_tag_list{{atom_tag, atom_tag}};

      std::array<size_t, PairLayer + 1> self_indices_pair;
      self_indices_pair[PairLayer] = pair_counter;
      auto && self_indices_pair_ = IndexConstArray(self_indices_pair.data());
      pair_cluster_indices.push_back(self_indices_pair_);
      auto && self_pair =
          ClusterRefKey<2, PairLayer>(self_atom_tag_list, self_indices_pair_);
      this->add_atom(self_pair);
      pair_counter++;

      for (auto pair : atom) {
        this->add_atom(pair);

        Eigen::Matrix<size_t, PairLayer + 1, 1> indices_pair;
        indices_pair(PairLayer) = pair_counter;
        pair_cluster_indices.push_back(indices_pair);
        pair_counter++;
      }
    }
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_CONTRIBUTION_HH_
