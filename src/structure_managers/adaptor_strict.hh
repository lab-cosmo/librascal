/**
 * file   adaptor_strict.hh
 *
 * @author Till Junge <till.junge@altermail.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   04 Jun 2018
 *
 * @brief implements an adaptor for structure_managers, filtering
 * the original manager so that only neighbours that are strictly
 * within r_cut are retained
 *
 * Copyright © 2018 Till Junge, Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_STRICT_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_STRICT_HH_

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "structure_managers/updateable_base.hh"
#include "rascal_utility.hh"

namespace rascal {
  /*
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorStrict;

  /*
   * specialisation of traits for strict adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorStrict<ManagerImplementation>> {
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::yes};
    constexpr static bool HasDistances{true};
    constexpr static bool HasDirectionVectors{true};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder};
    using LayerByOrder = typename LayerIncreaser<
        MaxOrder, typename ManagerImplementation::traits::LayerByOrder>::type;
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
  class AdaptorStrict
      : public StructureManager<AdaptorStrict<ManagerImplementation>>,
        public std::enable_shared_from_this<
            AdaptorStrict<ManagerImplementation>> {
   public:
    using Manager_t = AdaptorStrict<ManagerImplementation>;
    using Parent = StructureManager<Manager_t>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using traits = StructureManager_traits<AdaptorStrict>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;

    static_assert(traits::MaxOrder > 1,
                  "ManagerImlementation needs to handle pairs");

    //! Default constructor
    AdaptorStrict() = delete;

    /**
     * construct a strict neighbourhood list from a given manager. `cut-off`
     * specifies the strict cutoff radius. all clusters with distances above
     * this parameter will be skipped
     */
    AdaptorStrict(ImplementationPtr_t manager, double cutoff);

    AdaptorStrict(ImplementationPtr_t manager, std::tuple<double> tp)
        : AdaptorStrict(manager, std::get<0>(tp)) {}

    AdaptorStrict(ImplementationPtr_t manager, const Hypers_t & adaptor_hypers)
        : AdaptorStrict(manager,
                        adaptor_hypers.at("cutoff").template get<double>()) {}

    //! Copy constructor
    AdaptorStrict(const AdaptorStrict & other) = delete;

    //! Move constructor
    AdaptorStrict(AdaptorStrict && other) = default;

    //! Destructor
    virtual ~AdaptorStrict() = default;

    //! Copy assignment operator
    AdaptorStrict & operator=(const AdaptorStrict & other) = delete;

    //! Move assignment operator
    AdaptorStrict & operator=(AdaptorStrict && other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    inline void update_self();

    //! update the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments);

    //! returns the (strict) cutoff for the adaptor
    inline const double & get_cutoff() const { return this->cutoff; }

    //! returns the distance between atoms in a given pair
    template <size_t Order, size_t Layer>
    inline const double &
    get_distance(const ClusterRefKey<Order, Layer> & pair) const {
      return this->distance[pair];
    }

    template <size_t Order, size_t Layer>
    inline double & get_distance(const ClusterRefKey<Order, Layer> & pair) {
      return this->distance[pair];
    }

    //! returns the direction vector between atoms in a given pair
    template <size_t Order, size_t Layer>
    inline const Vector_ref
    get_direction_vector(const ClusterRefKey<Order, Layer> & pair) const {
      return this->dir_vec[pair];
    }

    template <size_t Order, size_t Layer>
    inline Vector_ref
    get_direction_vector(const ClusterRefKey<Order, Layer> & pair) {
      return this->dir_vec[pair];
    }

    inline size_t get_nb_clusters(int order) const {
      return this->atom_indices[order - 1].size();
    }

    inline size_t get_size() const { return this->manager->get_size(); }

    inline size_t get_size_with_ghosts() const {
      return this->get_nb_clusters(1);
    }

    inline Vector_ref get_position(const int & index) {
      return this->manager->get_position(index);
    }

    //! get atom_index of index-th neighbour of this cluster
    template <size_t Order, size_t Layer>
    inline int
    get_cluster_neighbour(const ClusterRefKey<Order, Layer> & cluster,
                          int index) const {
      static_assert(Order <= traits::MaxOrder - 1,
                    "this implementation only handles upto traits::MaxOrder");
      auto && offset = this->offsets[Order][cluster.get_cluster_index(Layer)];
      return this->atom_indices[Order][offset + index];
    }

    //! get atom_index of the index-th atom in manager
    inline int get_cluster_neighbour(const Parent & /*parent*/,
                                     size_t index) const {
      return this->atom_indices[0][index];
    }

    //! return atom type
    inline int & get_atom_type(const AtomRef_t & atom) {
      /**
       * careful, atom refers to our local index, for the manager, we need its
       * index:
       */
      auto && original_atom{this->atom_indices[0][atom.get_index()]};
      return this->manager->get_atom_type(original_atom);
    }

    //! return atom type
    inline const int & get_atom_type(const AtomRef_t & atom) const {
      // careful, atom refers to our local index, for the manager, we need its
      // index:
      auto && original_atom{this->atom_indices[0][atom.get_index()]};
      return this->manager->get_atom_type(original_atom);
    }

    //! Returns atom type given an atom index
    inline int & get_atom_type(const int & atom_id) {
      return this->manager->get_atom_type(atom_id);
    }

    //! Returns a constant atom type given an atom index
    inline const int & get_atom_type(const int & atom_id) const {
      auto && type{this->manager->get_atom_type(atom_id)};
      return type;
    }
    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template <size_t Order>
    inline size_t
    get_offset_impl(const std::array<size_t, Order> & counters) const {
      static_assert(Order < traits::MaxOrder,
                    "Calling this function with the wrong order cluster");
      return this->offsets[Order][counters.back()];
    }

    //! return the number of neighbours of a given atom
    template <size_t Order, size_t CallingLayer>
    inline size_t
    get_cluster_size(const ClusterRefKey<Order, CallingLayer> & cluster) const {
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

   protected:
    /**
     * main function during construction of a neighbourlist.
     * @param atom the atom to add to the list
     * @param Order select whether it is an i-atom (order=1), j-atom (order=2),
     * or ...
     */
    template <size_t Order>
    inline void add_atom(int atom_index) {
      static_assert(Order <= traits::MaxOrder,
                    "you can only add neighbours to the n-th degree defined by "
                    "MaxOrder of the underlying manager");

      // add new atom at this Order
      this->atom_indices[Order].push_back(atom_index);
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

    template <size_t Order>
    inline void
    add_atom(const typename ManagerImplementation::template ClusterRef<Order> &
                 cluster) {
      return this->template add_atom<Order - 1>(cluster.back());
    }

    ImplementationPtr_t manager;
    typename AdaptorStrict::template Property_t<double, 2> distance;
    typename AdaptorStrict::template Property_t<double, 2, 3> dir_vec;
    const double cutoff;

    /**
     * store atom indices per order,i.e.
     *   - atom_indices[0] lists all i-atoms
     *   - atom_indices[1] lists all j-atoms
     *   - atom_indices[2] lists all k-atoms
     *   - etc
     */
    std::array<std::vector<int>, traits::MaxOrder> atom_indices;
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

  namespace internal {
    /* ---------------------------------------------------------------------- */
    template <bool IsStrict, class ManagerImplementation>
    struct CutOffChecker {
      static bool check(const std::shared_ptr<ManagerImplementation> & manager,
                        double cutoff) {
        return cutoff < manager->get_cutoff();
      }
    };

    /* ---------------------------------------------------------------------- */
    template <class ManagerImplementation>
    struct CutOffChecker<false, ManagerImplementation> {
      static bool
      check(const std::shared_ptr<ManagerImplementation> & /*manager*/,
            double /*cutoff*/) {
        return true;
      }
    };

    /* ---------------------------------------------------------------------- */
    template <class ManagerImplementation>
    bool inline check_cutoff(
        const std::shared_ptr<ManagerImplementation> & manager, double cutoff) {
      constexpr bool IsStrict{(ManagerImplementation::traits::Strict ==
                               AdaptorTraits::Strict::yes)};
      return CutOffChecker<IsStrict, ManagerImplementation>::check(manager,
                                                                   cutoff);
    }
  }  // namespace internal

  /*--------------------------------------------------------------------------*/
  template <class ManagerImplementation>
  AdaptorStrict<ManagerImplementation>::AdaptorStrict(
      std::shared_ptr<ManagerImplementation> manager, double cutoff)
      : manager{std::move(manager)}, distance{*this}, dir_vec{*this},
        cutoff{cutoff}, atom_indices{}, nb_neigh{}, offsets{}

  {
    if (not internal::check_cutoff(this->manager, cutoff)) {
      throw std::runtime_error("underlying manager already has a smaller "
                               "cut off");
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class... Args>
  void AdaptorStrict<ManagerImplementation>::update(Args &&... arguments) {
    if (sizeof...(arguments) > 0) {
      this->set_update_status(false);
    }
    this->manager->update(std::forward<Args>(arguments)...);
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorStrict<ManagerImplementation>::update_self() {
    //! Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    //! initialise the neighbourlist
    for (size_t i{0}; i < traits::MaxOrder; ++i) {
      this->atom_indices[i].clear();
      this->nb_neigh[i].clear();
      this->offsets[i].clear();
    }
    this->nb_neigh[0].push_back(0);
    for (auto & vector : this->offsets) {
      vector.push_back(0);
    }

    //! initialise the distance storage
    this->distance.resize_to_zero();
    this->dir_vec.resize_to_zero();

    // fill the list, at least pairs are mandatory for this to work
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    size_t pair_counter{0};
    // depending on the underlying neighbourlist, the proxy `.with_ghosts()` is
    // either actually with ghosts, or only returns the number of centers.
    for (auto atom : this->manager.get()->with_ghosts()) {
      this->add_atom(atom);
      /**
       * Add new layer for atoms (see LayerByOrder for
       * possible optimisation).
       */

      constexpr auto AtomLayer{
          compute_cluster_layer<atom.order()>(typename traits::LayerByOrder{})};

      Eigen::Matrix<size_t, AtomLayer + 1, 1> indices;

      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer - 1);
      atom_cluster_indices.push_back(indices);
      double rc2{this->cutoff * this->cutoff};
      for (auto pair : atom) {
        constexpr auto PairLayer{compute_cluster_layer<pair.order()>(
            typename traits::LayerByOrder{})};

        auto vec_ij{pair.get_position() - atom.get_position()};
        double distance2{(vec_ij).squaredNorm()};

        if (distance2 <= rc2) {
          this->add_atom(pair);
          double distance{std::sqrt(distance2)};

          this->dir_vec.push_back((vec_ij.array() / distance).matrix());
          this->distance.push_back(distance);

          Eigen::Matrix<size_t, PairLayer + 1, 1> indices_pair;
          indices_pair.template head<PairLayer>() = pair.get_cluster_indices();
          indices_pair(PairLayer) = pair_counter;
          pair_cluster_indices.push_back(indices_pair);

          pair_counter++;
        }
      }
    }
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_STRICT_HH_
