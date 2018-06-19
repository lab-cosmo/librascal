/**
 * file   adaptor_increase_maxlevel.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   19 Jun 2018
 *
 * @brief implements an adaptor for neighbourhood_managers, which
 * creates a full and half neighbourlist if there is none and
 * triplets/quadruplets, etc. if existant.
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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


#ifndef ADAPTOR_MAXLEVEL_H
#define ADAPTOR_MAXLEVEL_H

namespace rascal {
  /**
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorMaxLevel;

  /**
   * specialisation of traits for increase MaxLevel adaptor
   */
  template <class ManagerImplementation>
  struct NeighbourhoodManager_traits<AdaptorMaxLevel<ManagerImplementation>> {

    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
      ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    constexpr static int MaxLevel{ManagerImplementation::traits::MaxLevel};
  };

  /**
   * Adaptor that increases the MaxLevel of an existing
   * NeighbourhoodManager. This means, if the manager does not have a
   * neighbourlist, it is created, if it exists, triplets,
   * quadruplets, etc. lists are created.
   */
  template <class ManagerImplementation>
  class AdaptorMaxLevel: public
  NeighbourhoodManagerBase<AdaptorMaxLevel<ManagerImplementation>>
  {
  public:
    using Parent = ManagerImplementation;
    using traits = NeighbourhoodManager_traits<AdaptorMaxLevel>;
    using AtomRef_t = typename Parent::AtomRef_t;
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template ClusterRef_t<Level, MaxLevel>;
    using PairRef_t = ClusterRef_t<2, traits::MaxLevel>;

    // TODO if MaxLevel can be == 1 -> neighbourlist need to be built.
    static_assert(traits::MaxLevel > 1,
                  "ManagerImlementation needs to handle pairs");

    //! Default constructor
    AdaptorMaxLevel() = delete;

    /**
     * construct a strict neighbourhood list from a given manager and cut-off radius
     */
    AdaptorMaxLevel(ManagerImplementation& manager, double cutoff);

    //! Copy constructor
    AdaptorMaxLevel(const AdaptorMaxLevel &other) = delete;

    //! Move constructor
    AdaptorMaxLevel(AdaptorMaxLevel &&other) = default;

    //! Destructor
    virtual ~AdaptorMaxLevel() = default;

    //! Copy assignment operator
    AdaptorMaxLevel& operator=(const AdaptorMaxLevel &other) = delete;

    //! Move assignment operator
    AdaptorMaxLevel& operator=(AdaptorMaxLevel &&other) = default;

    //! Update just the adaptor assuming the underlying manager was
    //! updated. this function invokes building either the neighbour
    //! list or to make triplets, quadruplets, etc. depending on the
    //! MaxLevel
    void update();

    //! Update the underlying manager as well as the adaptor
    template<class ... Args>
    void update(Args&&... arguments);

    inline double get_cutoff() const {return this->cutoff;}

    //! Get number of atoms, pairs, triplets, etc.
    inline size_t get_nb_clusters(int cluster_size) const {
      return this->atom_refs[cluster_size-1].size();
    }

    //! Return number of atoms
    inline size_t get_size() const {
      return this->get_nb_clusters(1);
    }

    inline Vector_ref get_position(const AtomRef_t & atom) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[0][atom.get_index()]};
      return this->manager.get_position(original_atom);
    }

    template<int Level, int MaxLevel>
    inline Vector_ref get_neighbour_position(const ClusterRef_t<Level,
					     MaxLevel>& cluster) {
      // Argument is now the same, but implementation
      throw std::runtime_error("should be adapted to Félix's new interface using the ClusterRef");
    }


    // return the global id of an atom
    inline size_t get_atom_id(const Parent& parent, int i_atom_id) const {
      return this->manager.get_atom_id(parent, i_atom_id);
    }

    // return the global id of an atom
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const {
      // careful, cluster refers to our local index, for the manager, we
      // need its index:
      auto original_cluster{cluster};
      auto & original_atoms = original_atoms.get_atoms();
      const auto &atoms = cluster.get_atoms();
      for (int i{0}; i < Level; i++) {
        original_atoms[i] = this->atom_refs[0][atoms[i].get_index()];
      }

      return this->manager.get_atom_id(original_cluster, j_atom_id);
    }


    //! Return atom type
    inline int & get_atom_type(const AtomRef_t& atom) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[0][atom.get_index()]};
      return this->manager.get_atom_type(original_atom);
    }

    //! Return atom type
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
    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const {
      return this->offsets[Level][cluster.get_index()];
    }

  protected:
    /**
     * main function during construction of a neighbourlist.  @param
     * atom the atom to add to the list @param level select whether it
     * is an i-atom (level=0), j-atom (level=1), or ...
     */
    template <int Level>
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

    template <int Level>
    inline void add_atom(typename ManagerImplementation::template
                         ClusterRef<Level, traits::MaxLevel> cluster) {
      return this->template add_atom <Level-1>(cluster.get_atoms().back());
    }

    template <int Level, bool IsDummy>
    struct HelperLoop;

    ManagerImplementation & manager;
    const double cutoff;

    /**
     * store atom refs per level if it exists,i.e.
     *   - atom_refs[0] lists all i-atoms
     *   - atom_refs[1] lists all j-atoms
     *   - atom_refs[2] lists all k-atoms
     *   - etc
     */
    // TODO: minimum is to have atoms from a neighbour list.
    std::array<std::vector<AtomRef_t>, traits::MaxLevel> atom_refs;
    // TODO these lists might need to be built, if a neighbourlist
    // does not exist.

    // /**
    //  * store the number of j-atoms for every i-atom (nb_neigh[1]), the
    //  * number of k-atoms for every j-atom (nb_neigh[2]), etc
    //  */
    std::array<std::vector<unsigned int>, traits::MaxLevel> nb_neigh;
    // /**
    //  * store the offsets from where the nb_neigh can be counted
    //  */
    std::array<std::vector<unsigned int>, traits::MaxLevel> offsets;
  private:
  };




  namespace internal {
    // TODO Add all neighbour list stuff in internal

    /* ---------------------------------------------------------------------- */



  }  // internal

  //----------------------------------------------------------------------------//
  // TODO include a distinction for the cutoff: with respect to the
  // i-atom only or with respect to the j,k,l etc. atom. I.e. the
  // cutoff goes with the Level.
  template <class ManagerImplementation>
  AdaptorMaxLevel<ManagerImplementation>::
  AdaptorMaxLevel(ManagerImplementation & manager, double cutoff):
    manager{manager},
    cutoff{cutoff},
    atom_refs{},
    nb_neigh{},
    offsets{}

  {
    if (traits::MaxLevel < 1) {
      throw std::runtime_error("No atoms in manager.");
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class ... Args>
  void AdaptorMaxLevel<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }

  /* ---------------------------------------------------------------------- */
  //! TODO

  void AdaptorMaxLevel::make_verlet_list() {
    // Make a verlet neighbour list according to Tadmor and Miller
    // 'Modeling Materials', algorithm 6.7, p 323, with modification
    // for periodicity. It results in a <code>strict</code>
    // neighbourlist.
    // This is only necessary, if the ManagerImplementation, with
    // which this adaptor is initializer does not have at least
    // already the atomic pairs.  Inluding ghost neighbours?  TODO:
    // add functionality for the shift vector

  }


  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxLevel<ManagerImplementation>::update() {
    // initialise the list if it does not exist
    // TODO
    // initialise the neighbourlist

    if (traits::MaxLevel == 1) {
      // Make full and half neighbour list (strict?)
    } else {
      // Make triplets/quadruplets/etc. based on existing
      // neighbourlist
      // Templated function?
    }

    // Initilize the list
    for (int i{0}; i < traits::MaxLevel; ++i) {
      this->atom_refs[i].clear();
      this->nb_neigh[i].resize(0);
      this->offsets[i].resize(0);
    }
    if (traits::MaxLevel > 1) {
      this->nb_neigh[0].push_back(0);
      for (auto & vector: this->offsets) {
	vector.push_back(0);
      }
    }
  }

  /* ---------------------------------------------------------------------- */



}  // rascal

#endif /* ADAPTOR_MAXLEVEL_H */
