/**
 * @file   rascal/structure_managers/adaptor_link_pairs.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   28 Jan 2020
 *
 * @brief implements an adaptor for structure_managers that allow for linking
 *        pairs ij with ji in a full neighborlist
 *
 * Copyright 2020 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_STRUCTURE_MANAGERS_LINK_PAIRS_HH_
#define SRC_RASCAL_STRUCTURE_MANAGERS_LINK_PAIRS_HH_

#include "rascal/structure_managers/property.hh"
#include "rascal/structure_managers/cluster_ref_key.hh"
#include "rascal/structure_managers/structure_manager.hh"
#include "rascal/structure_managers/updateable_base.hh"
#include "rascal/utils/utils.hh"

namespace rascal {
  /*
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class LinkPairs;

  /*
   * specialisation of traits for strict adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<LinkPairs<ManagerImplementation>> {
    using parent_traits = StructureManager_traits<ManagerImplementation>;
    constexpr static AdaptorTraits::Strict Strict{parent_traits::Strict};
    constexpr static bool HasDistances{parent_traits::HasDistances};
    constexpr static bool HasDirectionVectors{parent_traits::HasDirectionVectors};
    constexpr static bool HasSwapIJ{true};
    constexpr static bool HasCenterPair{parent_traits::HasCenterPair};
    constexpr static int Dim{parent_traits::Dim};
    constexpr static size_t MaxOrder{parent_traits::MaxOrder};
    constexpr static int StackLevel{parent_traits::StackLevel + 1};
    using LayerByOrder =
        typename LayerIncreaser<typename parent_traits::LayerByOrder>::type;
    using PreviousManager_t = ManagerImplementation;
    constexpr static AdaptorTraits::NeighbourListType NeighbourListType{
        parent_traits::NeighbourListType};
  };

  /**
   * Adaptor that guarantees that only neighbours within the cutoff are
   * present. A neighbor manager could include some wiggle room and list
   * clusters with distances above the specified cutoff, this adaptor makes it
   * possible to get a list with only the clusters that have distances strictly
   * below the cutoff. This is also useful to extract managers with different
   * levels of truncation from a single, loose manager.
   *
   * This interface should be implemented by all managers with the trait
   * AdaptorTraits::Strict::yes
   */
  template <class ManagerImplementation>
  class LinkPairs
      : public StructureManager<LinkPairs<ManagerImplementation>>,
        public std::enable_shared_from_this<
            LinkPairs<ManagerImplementation>> {
   public:
    using Manager_t = LinkPairs<ManagerImplementation>;
    using ManagerImplementation_t = ManagerImplementation;
    using Parent = StructureManager<Manager_t>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using traits = StructureManager_traits<LinkPairs>;
    using PreviousManager_t = typename traits::PreviousManager_t;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;
    using This = LinkPairs;
    using Distance_t = typename This::template Property_t<double, 2, 1>;
    using DirectionVector_t = typename This::template Property_t<double, 2, 3>;

    static_assert(traits::MaxOrder == 2,
                  "ManagerImlementation needs to handle pairs");
    static_assert(traits::NeighbourListType == AdaptorTraits::NeighbourListType::full,
                  "ManagerImlementation needs to be a full neighbor list.");
    static_assert(traits::HasDistances and traits::HasDirectionVectors,
          "ManagerImlementation needs to have distances and direction vectors");
    constexpr static auto AtomLayer{
        Manager_t::template cluster_layer_from_order<1>()};
    constexpr static auto PairLayer{
        Manager_t::template cluster_layer_from_order<2>()};

    //! Default constructor
    LinkPairs() = delete;

    /**
     */
    LinkPairs(ImplementationPtr_t manager);

    LinkPairs(ImplementationPtr_t manager, const Hypers_t & /*adaptor_hypers*/)
        : LinkPairs(manager) {}

    //! Copy constructor
    LinkPairs(const LinkPairs & other) = delete;

    //! Move constructor
    LinkPairs(LinkPairs && other) = default;

    //! Destructor
    virtual ~LinkPairs() = default;

    //! Copy assignment operator
    LinkPairs & operator=(const LinkPairs & other) = delete;

    //! Move assignment operator
    LinkPairs & operator=(LinkPairs && other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    void update_self();

    //! update the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments);

    //! returns the (strict) cutoff for the adaptor
    double get_cutoff() const { return this->manager->get_cutoff(); }

    size_t get_nb_clusters(int order) const {
      return this->manager->get_nb_clusters(order);
    }

    size_t get_size() const { return this->manager->get_size(); }

    Vector_ref get_position(int index) {
      return this->manager->get_position(index);
    }

    //! get atom_tag of index-th neighbour of this cluster
    template <size_t Order, size_t Layer>
    int get_neighbour_atom_tag(const ClusterRefKey<Order, Layer> & cluster,
                               int index) const {
      static_assert(Order < traits::MaxOrder,
                    "Calling this function with the wrong order cluster");
      return this->manager->get_neighbour_atom_tag(cluster, index);
    }

    //! get atom_tag of the index-th atom in manager
    int get_neighbour_atom_tag(const Parent & /*parent*/, size_t index) const {
      return this->center_tag_list[index];
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
      return this->manager->get_atom_type(atom);
    }

    //! Returns a constant atom type given an atom tag
    int get_atom_type(int atom_id) const {
      return this->manager->get_atom_type(atom_id);
    }
    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template <size_t Order>
    size_t get_offset_impl(const std::array<size_t, Order> & counters) const {
      return this->manager->get_offset_impl(counters);
    }

    template <size_t Order>
    int get_neighbour_atom_tag(const size_t neighbour_index) const {
      static_assert(Order < traits::MaxOrder,
                    "Calling this function with the wrong order cluster");
      return this->manager->get_neighbour_atom_tag(neighbour_index);
    }

    //! Returns the number of neighbours of a given atom at a given TargetOrder
    //! Returns the number of pairs of a given center
    template <size_t TargetOrder, size_t Order, size_t Layer>
    size_t get_cluster_size_impl(const ClusterRefKey<Order, Layer> & cluster) const {
      return this->manager->template get_cluster_size_impl<TargetOrder>(cluster);
    }

    //! Get the manager used to build the instance
    ImplementationPtr_t get_previous_manager_impl() {
      return this->manager->get_shared_ptr();
    }

    template <size_t Order, size_t Layer>
    size_t get_pair_ji_offset(const ClusterRefKey<Order, Layer> & pair) {
      static_assert(Order == 2, "This swap makes sense only with pairs.");
      auto && access_index = pair.get_cluster_index(PairLayer);
      return this->pairs_ji_offset.at(access_index);
    }

   protected:
    /**
     * Create the necessary data to associate the ij pair with the ji pair
     */
    void link_pairs();

    ImplementationPtr_t manager;

    /**
     *
     */
    std::vector<size_t> pairs_ji_offset;

    std::vector<int> center_tag_list;
  };

  /*--------------------------------------------------------------------------*/
  template <class ManagerImplementation>
  LinkPairs<ManagerImplementation>::LinkPairs(
      std::shared_ptr<ManagerImplementation> manager)
      : manager{std::move(manager)}, pairs_ji_offset{}, center_tag_list{} {}

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class... Args>
  void LinkPairs<ManagerImplementation>::update(Args &&... arguments) {
    this->manager->update(std::forward<Args>(arguments)...);
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void LinkPairs<ManagerImplementation>::link_pairs() {
    for (auto center : this->manager) {
      const int atom_i_tag{center.get_atom_tag()};
      for (auto pair_ij : center.pairs_with_self_pair()) {
        const double& dist_ij{this->manager->get_distance(pair_ij)};
        const auto dir_vec_ij = this->manager->get_direction_vector(pair_ij);
        const int atom_j_tag = pair_ij.get_internal_neighbour_atom_tag();
        const size_t atom_j_index = this->manager->get_atom_index(atom_j_tag);
        auto atom_j_it = this->manager->get_iterator_at(atom_j_index, 0);
        auto atom_j = *atom_j_it;
        size_t pair_offset{0};
        // loop over center j to find the ji pair corresponding the ij pair
        for (auto pair_ji : atom_j.pairs_with_self_pair()) {
          auto atom_i_prime = pair_ji.get_atom_j();
          const int atom_i_prime_tag = atom_i_prime.get_atom_tag();
          // check that atom_i_prime corresponds to either center i or one
          // of its periodic images
          if (atom_i_prime_tag == atom_i_tag) {
            const double& dist_ji{this->manager->get_distance(pair_ji)};
            // this should be a sufficient condition in most cases
            if (dist_ij - dist_ji < 1e-13) {
              const auto dir_vec_ji = this->manager->get_direction_vector(pair_ji);
              // direction vectors should be of oposite sign
              if (((dir_vec_ij + dir_vec_ji).array() < 1e-13).all()) {
                this->pairs_ji_offset.emplace_back(pair_offset);
                break;
              }
            }
          }
          ++pair_offset;
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void LinkPairs<ManagerImplementation>::update_self() {
    //! Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    this->pairs_ji_offset.clear();
    this->center_tag_list.clear();
    // fill the list, at least pairs are mandatory for this to work
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    size_t pair_counter{0};

    for (auto && atom : this->manager) {
      Eigen::Matrix<size_t, AtomLayer + 1, 1> indices;
      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer - 1);
      atom_cluster_indices.push_back(indices);
      this->center_tag_list.emplace_back(atom.get_atom_tag());
      for (auto pair : atom.pairs_with_self_pair()) {
        Eigen::Matrix<size_t, PairLayer + 1, 1> indices_pair;
        indices_pair.template head<PairLayer>() = pair.get_cluster_indices();
        indices_pair(PairLayer) = pair_counter;
        pair_cluster_indices.push_back(indices_pair);
        pair_counter++;
      }
    }

    for (auto && atom : this->manager->only_ghosts()) {
      Eigen::Matrix<size_t, AtomLayer + 1, 1> indices;
      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer - 1);
      atom_cluster_indices.push_back(indices);
      this->center_tag_list.emplace_back(atom.get_atom_tag());
    }

    this->link_pairs();
  }
}  // namespace rascal

#endif  // SRC_RASCAL_STRUCTURE_MANAGERS_LINK_PAIRS_HH_
