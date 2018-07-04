/**
 * file   adaptor_strict.hh
 *
 * @author Till Junge <till.junge@altermail.ch>
 *
 * @date   04 Jun 2018
 *
 * @brief implements an adaptor for neighbourhood_managers, filtering
 * the original manager so that only neighbours that are strictly
 * within r_cut are retained
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/property.hh"


#ifndef ADAPTOR_STRICT_H
#define ADAPTOR_STRICT_H

namespace rascal {
  /**
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorStrict;

  /**
   * specialisation of traits for strict adaptor
   */
  template <class ManagerImplementation>
  struct NeighbourhoodManager_traits<AdaptorStrict<ManagerImplementation>> {

    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::yes};
    constexpr static bool HasDistances{true};
    constexpr static bool HasDirectionVectors{ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    constexpr static size_t MaxLevel{ManagerImplementation::traits::MaxLevel};
    using Depth = DepthIncreaser<MaxLevel, ManagerImplementation::traits::Depth>::type;
  };

  /**
   * Adaptor that guarantees that only neighbours within the cutoff
   * are present. This interface should be implemented by all managers
   * with the trait AdaptorTraits::Strict::yes
   */
  template <class ManagerImplementation>
  class AdaptorStrict: public
  NeighbourhoodManagerBase<AdaptorStrict<ManagerImplementation>>
  {
  public:
    using Parent = NeighbourhoodManagerBase<AdaptorStrict<ManagerImplementation>>;
    using traits = NeighbourhoodManager_traits<AdaptorStrict>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Level>
    using ClusterRef_t = typename ManagerImplementation::template ClusterRef<Level>;
    using PairRef_t = ClusterRef_t<2>;

    static_assert(traits::MaxLevel > 1,
                  "ManagerImlementation needs to handle pairs");

    //! Default constructor
    AdaptorStrict() = delete;

    /**
     * construct a strict neighbourhood list from a given manager and cut-off radius
     */
    AdaptorStrict(ManagerImplementation& manager, double cut_off);

    //! Copy constructor
    AdaptorStrict(const AdaptorStrict &other) = delete;

    //! Move constructor
    AdaptorStrict(AdaptorStrict &&other) = default;

    //! Destructor
    virtual ~AdaptorStrict() = default;

    //! Copy assignment operator
    AdaptorStrict& operator=(const AdaptorStrict &other) = delete;

    //! Move assignment operator
    AdaptorStrict& operator=(AdaptorStrict &&other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    void update();

    //! update the underlying manager as well as the adaptor
    template<class ... Args>
    void update(Args&&... arguments);

    inline double get_cutoff() const {return this->cut_off;}

    inline const double & get_distance(const PairRef_t & pair) const {
      return this->distance[pair];
    }
    inline double & get_distance(const PairRef_t & pair) {
      return this->distance[pair];
    }

    inline size_t get_nb_clusters(int cluster_size) const {
      return this->atom_refs[cluster_size-1].size();
    }

    inline size_t get_size() const {
      return this->get_nb_clusters(1);
    }

    inline Vector_ref get_position(const AtomRef_t & atom) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[0][atom.get_index()]};
      return this->manager.get_position(original_atom);
    }

    template<size_t Level>
    inline Vector_ref get_neighbour_position(const ClusterRef_t<Level>&
                                             /*cluster*/) {
      static_assert(Level > 1,
                    "Only possible for Level > 1.");
      static_assert(Level <= traits::MaxLevel,
                    "this implementation should only work up to MaxLevel.");
      // Argument is now the same, but implementation
      throw std::runtime_error("should be adapted to Félix's new interface using the ClusterRef");

    }

    // return the global id of an atom
    inline size_t get_atom_id(const Parent& /*parent*/, int i_atom_id) const {
      return this->manager.get_atom_id(this->manager, i_atom_id);
    }

    // return the global id of an atom
    template<size_t Level>
    inline size_t get_atom_id(const ClusterRefBase<Level>& cluster,
                              int j_atom_id) const {
      static_assert(Level <= traits::MaxLevel-1,
                    "this implementation only handles upto traits::MaxLevel");
      return this->manager.get_atom_id(cluster, j_atom_id);
    }


    //! return atom type
    inline int & get_atom_type(const AtomRef_t& atom) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[0][atom.get_index()]};
      return this->manager.get_atom_type(original_atom);
    }

    //! return atom type
    inline const int & get_atom_type(const AtomRef_t& atom) const {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[0][atom.get_index()]};
      return this->manager.get_atom_type(original_atom);
    }

    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<size_t Level>
    inline int get_offset_impl(const ClusterRef_t<Level>& cluster) const {
      return this->offsets[Level][cluster.get_index()];
    }

    // return the number of neighbours of a given atom
    template<size_t Level>
    inline size_t get_cluster_size(const ClusterRefBase<Level>& cluster) const {
      static_assert(Level <= traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      return this->nb_neigh[Level][cluster.back()];
    }

  protected:
    /**
     * main function during construction of a neighbourlist.  @param
     * atom the atom to add to the list @param level select whether it
     * is an i-atom (level=0), j-atom (level=1), or ...
     */
    template <size_t Level>
    inline void add_atom(typename ManagerImplementation::AtomRef_t atom) {
      static_assert(Level <= traits::MaxLevel,
                    "you can only add neighbours to the n-th degree defined by "
                    "MaxLevel of the underlying manager");

      // add new atom at this Level
      this->atom_refs[Level].push_back(atom);
      // count that this atom is a new neighbour
      this->nb_neigh[Level].back()++;
      this->offsets[Level].back()++;

      for (int i{Level+1}; i < traits::MaxLevel; ++i) {
        // make sure that this atom starts with zero lower-Level neighbours
        this->nb_neigh[i].push_back(0);
        // update the offsets
        this->offsets[i].push_back(this->offsets[i].back() +
                                   this->nb_neigh[i-1].back());
      }
    }

    template <size_t Level>
    inline void add_atom(typename ManagerImplementation::template
                         ClusterRef<Level> cluster) {
      return this->template add_atom <Level-1>(cluster.back());
    }

    template <size_t Level, bool IsDummy>
    struct HelperLoop;

    ManagerImplementation & manager;
    Property<AdaptorStrict, double, 2> distance;
    const double cut_off;

    /**
     * store atom refs per level,i.e.
     *   - atom_refs[0] lists all i-atoms
     *   - atom_refs[1] lists all j-atoms
     *   - atom_refs[2] lists all k-atoms
     *   - etc
     */
    std::array<std::vector<AtomRef_t>, traits::MaxLevel> atom_refs;
    /**
     * store the number of j-atoms for every i-atom (nb_neigh[1]), the
     * number of k-atoms for every j-atom (nb_neigh[2]), etc
     */
    std::array<std::vector<unsigned int>, traits::MaxLevel> nb_neigh;
    /**
     * store the offsets from where the nb_neigh can be counted
     */
    std::array<std::vector<unsigned int>, traits::MaxLevel> offsets;
  private:
  };




  namespace internal {

    /* ---------------------------------------------------------------------- */
    template<bool IsStrict, class ManagerImplementation>
    struct CutOffChecker {
      static bool check(const ManagerImplementation & manager,
                        double cut_off) {
        return cut_off < manager.get_cutoff();
      }
    };

    /* ---------------------------------------------------------------------- */
    template<class ManagerImplementation>
    struct CutOffChecker<false, ManagerImplementation> {
      static bool check(const ManagerImplementation & /*manager*/,
                        double /*cut_off*/) {
        return true;
      }
    };

    /* ---------------------------------------------------------------------- */
    template <class ManagerImplementation>
    bool inline check_cut_off(const ManagerImplementation & manager,
                              double cut_off) {
      constexpr bool IsStrict{(ManagerImplementation::traits::Strict ==
                               AdaptorTraits::Strict::yes)};
      return CutOffChecker<IsStrict, ManagerImplementation>::
        check(manager, cut_off);
    }


  }  // internal

  //----------------------------------------------------------------------------//
  template <class ManagerImplementation>
  AdaptorStrict<ManagerImplementation>::
  AdaptorStrict(ManagerImplementation & manager, double cut_off):
    manager{manager},
    distance{*this},
    cut_off{cut_off},
    atom_refs{},
    nb_neigh{},
    offsets{}

  {
    if (not internal::check_cut_off(manager, cut_off)) {
      throw std::runtime_error("underlying manager already has a smaller cut_off");
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class ... Args>
  void AdaptorStrict<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }


  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <size_t Level, bool IsDummy>
  struct AdaptorStrict<ManagerImplementation>::HelperLoop {
    static constexpr size_t MaxLevel{ManagerImplementation::traits::MaxLevel};
    using ClusterRef_t = typename ManagerImplementation::template
      ClusterRef<Level>;

    using NextLevelLoop = HelperLoop<Level+1,
                                     (Level+1 == MaxLevel)>;

    static void loop(ClusterRef_t & cluster, AdaptorStrict& manager) {
      for (auto next_cluster: cluster) {
        manager.add_atom(next_cluster);
        NextLevelLoop::loop(next_cluster, manager);
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  /**
   * End of recursion
   */
  template <class ManagerImplementation>
  template <size_t Level>
  struct AdaptorStrict<ManagerImplementation>::HelperLoop<Level, true> {
    static constexpr size_t MaxLevel{ManagerImplementation::traits::MaxLevel};
    using ClusterRef_t = typename ManagerImplementation::template
      ClusterRef<Level>;
    static void loop(ClusterRef_t & /*cluster*/,
                     AdaptorStrict<ManagerImplementation>& /*manager*/) {
      // do nothing
    }
  };


  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorStrict<ManagerImplementation>::update() {
    // initialise the neighbourlist
    for (int i{0}; i < traits::MaxLevel; ++i) {
      this->atom_refs[i].clear();
      this->nb_neigh[i].resize(0);
      this->offsets[i].resize(0);
    }
    this->nb_neigh[0].push_back(0);
    for (auto & vector: this->offsets) {
      vector.push_back(0);
    }

    // initialise the distance storage
    this->distance.resize_to_zero();

    // fill the list
    for (auto atom: this->manager) {
      this->add_atom(atom);
      for (auto pair: atom) {
        double distance{(atom.get_position()-
                         pair.get_position()).norm()};
        if (distance <= this->cut_off) {
          this->add_atom(pair);
          this->distance.push_back(distance);
        }
        using HelperLoop = HelperLoop<2, 2 >= traits::MaxLevel>;
        HelperLoop::loop(pair, *this);
      }
    }
  }

  /* ---------------------------------------------------------------------- */



}  // rascal

#endif /* ADAPTOR_STRICT_H */
