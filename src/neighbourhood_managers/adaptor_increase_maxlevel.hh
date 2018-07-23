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
#include "rascal_utility.hh"

#include <typeinfo>
#include <set>

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
      DepthExtender<MaxLevel,
                    typename
                    ManagerImplementation::traits::DepthByDimension>::type;
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

    using Parent =
      NeighbourhoodManagerBase<AdaptorMaxLevel<ManagerImplementation>>;
    using traits = NeighbourhoodManager_traits<AdaptorMaxLevel>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Level>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Level>;
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

    template<size_t Level>
    inline size_t get_offset_impl(const std::array<size_t, Level>
				  & counters) const;

    //! Get number of MaxLevel+1 tuplets
    inline size_t get_nb_clusters() const {
      return this->atom_refs.size();
    }

    //! Return number of atoms
    inline size_t get_size() const {
      // return this->get_nb_clusters(1); TODO: this has to return the
      // original number of atoms from the underlying manager
      return this->manager.get_size();
    }

    inline Vector_ref get_position(const size_t & atom_index) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto && original_atom{this->atom_refs[atom_index]};
      return this->manager.get_position(original_atom);
    }

    inline Vector_ref get_position(const AtomRef_t & atom) {
      // careful, atom refers to our local index, for the manager, we
      // need its index:
      auto atom_index{atom.get_index()};
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
    inline int get_cluster_neighbour(const Parent& /*parent*/,
				     size_t index) const {
      return this->manager.get_cluster_neighbour(this->manager, index);
    }

    template<size_t Level, size_t Depth>
    inline int get_cluster_neighbour(const ClusterRefBase<Level, Depth>
				     & cluster,
				     size_t index) const {
      static_assert(Level < traits::MaxLevel,
                    "this implementation only handles upto traits::MaxLevel");
      if (Level < traits::MaxLevel-1) {
	return this->manager.get_cluster_neighbour(cluster, index);
      } else {
	auto && offset = this->offsets[cluster.get_cluster_index(Depth)];
	return this->neighbours[offset + index];
      }
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
    inline size_t get_cluster_size(const ClusterRefBase<Level, Depth>
                                   & cluster) const {
      static_assert(Level < traits::MaxLevel,
                    "this implementation handles only the respective MaxLevel");
      // TODO: Should this not be a test for Depth?
      if (Level < traits::MaxLevel-1) {
	return this->manager.get_cluster_size(cluster);
      } else {
	return nb_neigh[cluster.back()] ;
      }
      // return this->manager.get_cluster_size(cluster); TODO: this should be
      // self-referring if Level=MaxLevel, otherwise underlying manager
      // TODO: conditional_t ??
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
     * atom the atom to add to the list. since the MaxLevel is
     * increased by one in this adaptor, the Level=MaxLevel
     */

    inline void add_atom(const int atom_index) {
      // static_assert(Level <= traits::MaxLevel,
      //               "you can only add neighbours to the n-th degree defined by "
      //               "MaxLevel of the underlying manager");

      // add new atom at this Level
      this->atom_indices.push_back(atom_index);
      // count that this atom is a new neighbour
      this->nb_neigh.back()++;
      this->offsets.back()++;

      // make sure that this atom starts with zero lower-Level neighbours
      this->nb_neigh.push_back(0);
      // update the offsets
      this->offsets.push_back(this->offsets.back() +
                              this->nb_neigh.back());
    }

    template <size_t Level>
    inline void add_atom(const typename ManagerImplementation::template
                         ClusterRef<Level> & cluster) {
      std::cout << "add_atom (cluster) "
		<< cluster.back()
		<< " Level " << Level
		<< std::endl;
      return this->add_atom(cluster.back());
    }

    template <size_t Level>
    inline void add_atom(typename ManagerImplementation::template
                         ClusterRef<Level> & cluster) {
      std::cout << "            adding_atom_level_up" << std::endl;
      return this->template add_atom(cluster.back());
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

    // not necessary any more
    // // stores AtomRefs to of neighbours for traits::MaxLevel-1-*plets
    // std::vector<AtomRef_t> atom_refs{};

    // Stores atom indices
    std::vector<size_t> atom_indices{}; //akin to ilist[]

    // stores the number of neighbours of traits::MaxLevel-1-*plets
    std::vector<size_t> nb_neigh{};

    // stores all neighbours of traits::MaxLevel-1-*plets
    std::vector<size_t> neighbours{};

    // stores the offsets of traits::MaxLevel-1-*plets for accessing
    // <code>neighbours</code>
    std::vector<size_t> offsets{};

    size_t cluster_counter{0};
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
    atom_indices{},
    nb_neigh{0},
    offsets{0}

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
      // size_t cluster_counter{0};

      // do nothing if MaxLevel is not reached, except call the next level
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
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Level>;
    using ClusterRefZero_t
    =   typename ManagerImplementation::template ClusterRef<1>;
    using traits = typename AdaptorMaxLevel<ManagerImplementation>::traits;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;

    static void loop (ClusterRef_t & cluster,
		      AdaptorMaxLevel<ManagerImplementation> & manager) {

      // TODO: Understand why manager.get_cluster_indices(), compare to strict?!
      auto & next_cluster_indices{std::get<Level>(manager.get_cluster_indices())};

      // get new depth of cluster indices
      constexpr auto NextClusterDepth{
        compute_cluster_depth<cluster.level()>
          (typename traits::DepthByDimension{})};

      Eigen::Matrix<size_t, NextClusterDepth+1, 1> indices_cluster;

      std::cout << " ------At old MaxLevel, adding new layer----- "
                << OldMaxLevel << std::endl;

      // get all i_atom 'names' to find neighbours for
      auto i_atoms = cluster.get_atom_indices();
      std::cout << "no of atoms in cluster " << i_atoms.size() << std::endl;
      std::cout << "atom indices " << i_atoms[0] << " " << i_atoms[1] << std::endl;
      std::cout << "-------------------->" << std::endl;
      std::cout << "cluster.get_cluster_index(depth) "
                << cluster.get_cluster_index(0) << std::endl;
      std::cout << " Level " << cluster.level() << std::endl;

      ClusterRefZero_t a{cluster.front(), cluster.get_manager()};
      std::cout << "get_cluster_size " << cluster.get_manager().cluster_size(a) << std::endl;


      // // set for storing atom indices which have to be iterated for neighbours
      // std::set<int> current_i_atoms;

      // for (auto idx : i_atoms) {
      //   current_i_atoms.insert(idx);
      // };

      // translate i atom 'names' into cluster_indices for accessing neighbour
      // offsets
      std::set<size_t> current_j_atoms;


      /**
       * The following is the naive approach to what is necessary here. But
       * instead of looping over all of the atoms in the manager, we want to
       * build an iterator, which gives all neighbours of an atom
       */
      // for (auto i_atom: i_atoms) {
      //   for (auto atom: this->manager){ // clusterref
      //     if (atom.get_atom_index() == i_atom) {
      //       for (auto pair: atom) {
      //         do stuff;
      //       }
      //     }
      //   }
      // }

      // Get the corresponding neighbours of all current i atoms, i.e. build
      // pairs with the collected i_atoms




      // i refers to the local atom for which the pairs j are added as
      // MaxLevel+1
      // auto atom_i = cluster.back();
      // std::cout << "Atom i " << atom_i << std::endl;
      // std::cout << manager.get_size() << std::endl;
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

      std::cout << "<--------------------" << std::endl;
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
      // Level 1, Level variable is at 0, atoms, index 0
      // std::cout << "## for atoms in manager ##" << std::endl;
      using AddLevelLoop = AddLevelLoop<atom.level(),
                                        atom.level() == traits::MaxLevel-1>;
      AddLevelLoop::loop(atom, *this);
    }


    /**
     * To make triplets
     * Iteration for i-j adding i-k and j-k pairs
     * Iteration for i-j-k adding i-j and j-k and k-l pairs
     * etc.
     */

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

  /* ---------------------------------------------------------------------- */
  template<class ManagerImplementation>
  template<size_t Level>
  inline size_t AdaptorMaxLevel<ManagerImplementation>::
  get_offset_impl(const std::array<size_t, Level> & counters) const {

    static_assert(Level < traits::MaxLevel,
                  "this implementation handles only the respective MaxLevel");
    /**
     * Level accessor: 0 - atoms
     *                 1 - pairs
     *                 2 - triplets
     *                 etc.
     */

    if (Level == traits::MaxLevel-1) {
      /**
       * Reinterpret counters as a smaller array to call parent offset
       * multiplet. This can then be used to access the actual offset here.
       */
      std::cout << "IncMax if Level " << Level << std::endl;
      std::cout << "counters "
                << counters[0] << " "
                << counters[1] << " "
                << "size " << counters.size()
                << std::endl;
      /**
       * reinterpret_cast cuts off the last index in the cluster to access lower
       * level manager get_offset_impl
       */
      auto call_counters =
        reinterpret_cast<const std::array<size_t, Level-1> &>(counters);
      std::cout << "call_counters "
                << call_counters[0] << " "
                << "size " << call_counters.size()
                << std::endl;

      auto i{this->manager.get_offset_impl(call_counters)};

      auto j{counters.back()};
      auto main_offset{this->offsets[i]};
      return main_offset + j;
    } else {
      /**
       * If not accessible at this level, call lower Level offsets from lower
       * level manager(s).
       */
      std::cout << "IncMax else Level " << Level << std::endl;
      // return this->manager.get_offset_impl(counters);
    }
  }
}  // rascal

#endif /* ADAPTOR_MAXLEVEL_H */
