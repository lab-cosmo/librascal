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

#include <typeinfo>

#ifndef ADAPTOR_MAXLEVEL_H
#define ADAPTOR_MAXLEVEL_H

namespace rascal {
  /**
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorMaxLevel;

  /**
   * specialisation of traits for increase <code>MaxLevel</code> adaptor
   */
  template <class ManagerImplementation>
  struct NeighbourhoodManager_traits<AdaptorMaxLevel<ManagerImplementation>> {

    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
      ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    // New MaxLevel upon construction!
    constexpr static size_t MaxLevel{ManagerImplementation::traits::MaxLevel+1};
    // New Depth
    // TODO: Is this the correct way to initialize the increased depth?
    using DepthByDimension = typename
      DepthIncreaser<MaxLevel, typename ManagerImplementation::traits::DepthByDimension>::type;
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
    using Base = NeighbourhoodManagerBase<AdaptorMaxLevel<ManagerImplementation>>;

    using Parent = NeighbourhoodManagerBase<AdaptorMaxLevel<ManagerImplementation>>;
    using traits = NeighbourhoodManager_traits<AdaptorMaxLevel>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Level>
    using ClusterRef_t = typename ManagerImplementation::template ClusterRef<Level>;
    //using PairRef_t = ClusterRef_t<2, traits::MaxLevel>;

    // TODO if MaxLevel can be == 1 -> neighbourlist need to be built.
    static_assert(traits::MaxLevel > 1,
                  "ManagerImplementation needs to have an atom list.");

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

    //! Get number of MaxLevel+1 tuplets
    inline size_t get_nb_clusters(int cluster_size) const {
      return this->atom_refs.size();
    }

    //! Return number of atoms
    inline size_t get_size() const {
      // return this->get_nb_clusters(1); TODO: this has to return the
      // original number of atoms from the underlying manager
      return this->manager.get_size();
    }

    inline Vector_ref get_position(const int & atom_index) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[atom_index]};
      return this->manager.get_position(original_atom);
    }

    template<size_t Level, size_t Depth>
    inline Vector_ref get_neighbour_position(const ClusterRefBase<Level, Depth>&
                                             /*cluster*/) {
      static_assert(Level > 1,
                    "Only possible for Level > 1.");
      static_assert(Level <= traits::MaxLevel,
                    "this implementation should only work up to MaxLevel.");
      // Argument is now the same, but implementation
      throw std::runtime_error("should be adapted to Félix's "
                               "new interface using the ClusterRef");
    }

    // return the global id of an atom
    inline size_t get_atom_id(const Parent& /*parent*/, int i_atom_id) const {
      return this->manager.get_atom_id(this->manager, i_atom_id);
    }

    // return the global id of an atom
    template<size_t Level>
    inline size_t get_atom_id(const ClusterRef_t<Level>& cluster,
                              int j_atom_id) const {
      static_assert(Level < traits::MaxLevel,
                    "this implementation only handles upto traits::MaxLevel");
      return this->manager.get_atom_id(cluster, j_atom_id);
    }


    //! Return atom type
    inline int & get_atom_type(const AtomRef_t& atom) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[atom.get_index()]};
      return this->manager.get_atom_type(original_atom);
    }

    //! Return atom type
    inline const int & get_atom_type(const AtomRef_t& atom) const {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[atom.get_index()]};
      return this->manager.get_atom_type(original_atom);
    }

    template<size_t Level, size_t Depth>
    inline size_t get_cluster_size(const ClusterRefBase<Level, Depth>& cluster) const {
      static_assert(Level < traits::MaxLevel,
                    "this implementation only handles atoms and pairs");
      return this->manager.get_cluster_size(cluster);
    }

    // TODO: get a function answer highest Level cluster size = MaxLevel here
    // TODO2: or probably manager_base? 
    // template<size_t Level, size_t Depth>
    // inline size_t get_cluster_size(const ClusterRefBase<Level, Depth>& cluster) const {
    //   static_assert(Level < traits::MaxLevel,
    //                 "this implementation only handles atoms and pairs");
    //   return this->manager.get_cluster_size(cluster);
    // }
    
  protected:
    /**
     * main function during construction of a neighbourlist.  @param
     * atom the atom to add to the list @param level select whether it
     * is an i-atom (level=0), j-atom (level=1), or ...
     */
    // TODO: change add atom from AtomRef_t to size_t, because all access is with atom index?!
    template <size_t Level>
    inline void add_atom(typename ManagerImplementation::AtomRef_t atom) {
      static_assert(Level < traits::MaxLevel,
                    "you can only add neighbours to the n-th degree defined by "
                    "MaxLevel of the underlying manager");

      auto p = &atom;
      std::cout << p << std::endl;
      // add new atom at this Level
      // this->atom_refs.push_back(atom);
      // count that this atom is a new neighbour
      // if (Level == (traits::MaxLevel-1)) {
      // 	std::cout << "incr nb_neigh " << nb_neigh[Level].back() << std::endl;
      // }
      this->nb_neigh.back()++;
      this->offsets.back()++;

      std::cout << "inside add atom: nb_neigh[Level].back() "
		<< " Level " << Level << ", back "
		<< nb_neigh.back()
		<< std::endl;
    }

    template <size_t Level>
    inline void add_atom(typename ManagerImplementation::template
                         ClusterRef<Level> cluster) {
      std::cout << "add_atom (cluster) "
		<< cluster.back()
		<< " Level " << Level
		<< std::endl;
      return this->template add_atom <Level-1>(cluster.back());
    }

    template <size_t Level>
    inline void add_atom_level_up(typename ManagerImplementation::template
				  ClusterRef<Level> cluster) {

      std::cout << "            adding_atom_level_up" << std::endl;
      return this->template add_atom <Level>(cluster);
    }

    // Make a half neighbour list, by construction only Level=1 is supplied.
    void make_half_neighbour_list();

    // Make full list
    void make_full_neighbour_list();

    // Increase whatever level is present
    void increase_maxlevel();

    // Construction to get the correct ClusterRef



    ManagerImplementation & manager;
    const double cutoff;

    template<size_t Level, bool IsDummy> struct AddLevelLoop;

    // stores AtomRefs to of neighbours for traits::MaxLevel-1-*plets
    std::vector<AtomRef_t> atom_refs;

    // stores the number of neighbours of traits::MaxLevel-1-*plets
    std::vector<unsigned int> nb_neigh{};

    // stores all neighbours of traits::MaxLevel-1-*plets
    std::vector<unsigned int> neighbours{};

    // stores the offsets of traits::MaxLevel-1-*plets for accessing
    // <code>neighbours</code>
    std::vector<unsigned int> offsets{};
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
  template <class ManagerImplementation>
  template <size_t Level, bool IsDummy>
  struct AdaptorMaxLevel<ManagerImplementation>::AddLevelLoop {

    static constexpr int OldMaxLevel{ManagerImplementation::traits::MaxLevel};
    using ClusterRef_t = typename ManagerImplementation::template
      ClusterRef<Level>;

    using NextLevelLoop = AddLevelLoop<Level+1,
				       (Level+1 == OldMaxLevel)>;
    static void loop(ClusterRef_t & cluster,
		     AdaptorMaxLevel<ManagerImplementation>& manager) {
      // do nothing if MaxLevel is not reached, except call the next
      // level
      for (auto next_cluster : cluster) {
	NextLevelLoop::loop(next_cluster, manager);
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  // At desired MaxLevel (plus one), here is where the magic happens
  // and the neighbours of the same level are added as the Level+1.
  // add check for non half neighbour list
  template <class ManagerImplementation>
  template <size_t Level>
  struct AdaptorMaxLevel<ManagerImplementation>::AddLevelLoop<Level, true> {
    static constexpr int OldMaxLevel{ManagerImplementation::traits::MaxLevel};
    using ClusterRef_t = typename ManagerImplementation::template
      ClusterRef<Level>;

    static void loop (ClusterRef_t & cluster,
		      AdaptorMaxLevel<ManagerImplementation>& manager) {
      std::cout << " ------At old MaxLevel, adding new layer-----" << std::endl;

      // i refers to the local atom for which the pairs j are added as
      // MaxLevel+1
      auto atom_i = cluster.back();
      std::cout << "Atom i " << atom_i << std::endl;
      std::cout << manager.get_size() << std::endl;
      // auto index {cluster.get_global_index()};
      // for (int i{0};  i  < this->nb_neigh[index]; ++i) {
      // 	auto jndex{this->neighbours[this->offsets[index]+i]}
      // 	auto this->
      // }



      // std::cout << "atom_i " << atom_i << std::endl;

      // auto nneigh_i = ;
      // std::cout << "neigh " << neigh << std::endl;

      // for (auto next_cluster : cluster) {
      // 	std::cout << ">>next_cluster ----- @MaxLevel " << Level<< " index "
      // 		  << next_cluster.back() << std::endl;
      // 	manager.add_atom_level_up(next_cluster);
      // }
    }
    // end
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxLevel<ManagerImplementation>::make_half_neighbour_list() {
    // Make a half neighbour list (not quite Verlet, because of the
    // missing skin) according to Tadmor and Miller 'Modeling
    // Materials', algorithm 6.7, p 323, (needs modification for
    // periodicity).

    // It results in a <code>strict</code> neighbourlist. It is
    // only a half-neighbour list.  The most obvious difference is
    // that no 'skin' is used in conjunction with the cutoff.  This is
    // only necessary, if the ManagerImplementation, with which this
    // adaptor is initialized does not have at least already atomic
    // pairs. Inluding ghost neighbours?

    // TODO: add functionality for the shift vector??!

    // The zeroth level does not have neighbours
    this->nb_neigh.push_back(0);

    unsigned int nneigh_off{0};

    for (auto it=this->manager.begin(); it!=--this->manager.end(); ++it){
      // Add atom at this level this is just the standard list.
      auto atom_i = *it;
      this->add_atom(atom_i);

      auto jt = it;
      ++jt; // go to next atom in manager (no self-neighbour)
      for (;jt!=manager.end(); ++jt){
      	auto atom_j = *jt;
      	double distance{(atom_i.get_position() -
			 atom_j.get_position()).norm()};
      	if (distance <= this->cutoff) {
      	  // Store atom_j in neighbourlist of atom_i
      	  this->nb_neigh.back()++;
      	  nneigh_off += 1;
      	}
      }
      this->offsets.push_back(nneigh_off);
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxLevel<ManagerImplementation>::make_full_neighbour_list() {
    // Make a full neighbourlist, whithout fancy linked list or
    // cell. Also missing are periodic boundary conditions.

    // The zeroth level does not have neighbours
    this->nb_neigh.push_back(0);

    unsigned int nneigh_off{0};

    for (auto atom_i : this->manager){
      this->add_atom(atom_i);
      for (auto atom_j : this->manager){
      	if(&atom_i != &atom_j) { // avoid self-neighbouring
	  // TODO: this .get_position() function will always give the
	  // real position, never the shifted one.
      	  double distance{(atom_i.get_position() -
			   atom_j.get_position()).norm()};
      	  if (distance <= this->cutoff) {
      	    // Store atom_j in neighbourlist of atom_i
	    // this->atom_refs[1].push_back(atom_j.back());
	    this->add_atom_level_up(atom_j);
      	    this->nb_neigh.back()++;
      	    nneigh_off += 1;
      	  }
      	}
      }
      this->offsets.push_back(nneigh_off);
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxLevel<ManagerImplementation>::increase_maxlevel() {
    // Depending on the existing Level, this function increases the
    // MaxLevel by one by adding all neighbours of the last atom at
    // the end of the chain as new end of a tuple.
    // This results in each triplet is only existing once.

    // Attention, <code>traits::MaxLevel</code> is already increased
    // in the traits upon construction, therefore the MaxLevel needs
    // to be larger than 2 (i.e. a NeighbourhoodManager with a
    // pairlist is present to call this function here.)
    static_assert(traits::MaxLevel > 2, "No neighbourlist present.");

    for (auto atom : this->manager) {
      // Level 1, atoms, index 0
      std::cout << "## for atoms in manager ##" << std::endl;

      using AddLevelLoop = AddLevelLoop<1, 1 == traits::MaxLevel-1>;
      AddLevelLoop::loop(atom, *this);
    }


    // Iteration for i-j adding j-k pairs
    // Iteration for i-j-k addint k-l pairs

    // triplets:
    // atom_refs[new_max_level].push_back(j,k atoms)
    // quadruplets
    // atom_refs[new_max_level].push_back(j,k,l atoms)
    // The new level atom_refs need to get the correct ClusterRef to
    // be acessible as triplet/quadruplet/etc.
    // The existing MaxLevel is the iteration depth

  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxLevel<ManagerImplementation>::update() {
    // initialise the list if it does not exist
    // TODO
    // initialise the neighbourlist

    // -1 because the traits' MaxLevel is already increased
    if (traits::MaxLevel-1 == 1) {
      // Make half neighbour list (strict?)
      // initialise the neighbourlist
      for (int i{0}; i < traits::MaxLevel; ++i) {
	// this->atom_refs.clear();
	this->nb_neigh.resize(0);
	this->offsets.resize(0);
      }
      this->make_half_neighbour_list();
      // this->make_full_neighbour_list(); // no frills, full neighbourlist
    } else {
      // Make triplets/quadruplets/etc. based on existing
      // neighbourlist
      // Templated function?
      this->increase_maxlevel();
    }

  }

}  // rascal

#endif /* ADAPTOR_MAXLEVEL_H */
