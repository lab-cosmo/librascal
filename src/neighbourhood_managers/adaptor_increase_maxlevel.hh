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
   * Forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorMaxOrder;

  /**
   * Specialisation of traits for increase <code>MaxOrder</code> adaptor
   */
  template <class ManagerImplementation>
  struct NeighbourhoodManager_traits<AdaptorMaxOrder<ManagerImplementation>> {

    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
      ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    //! New MaxOrder upon construction!
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder+1};
    //! New Layer
    //! TODO: Is this the correct way to initialize the increased depth?
    using LayerByDimension = typename
      LayerExtender<MaxOrder,
                    typename
                    ManagerImplementation::traits::LayerByDimension>::type;
  };

  /**
   * Adaptor that increases the MaxOrder of an existing
   * NeighbourhoodManager. This means, if the manager does not have a
   * neighbourlist, it is created, if it exists, triplets, quadruplets,
   * etc. lists are created.
   */
  template <class ManagerImplementation>
  class AdaptorMaxOrder: public
  NeighbourhoodManagerBase<AdaptorMaxOrder<ManagerImplementation>>
  {
  public:
    using Base = NeighbourhoodManagerBase<AdaptorMaxOrder<ManagerImplementation>>;

    using Parent =
      NeighbourhoodManagerBase<AdaptorMaxOrder<ManagerImplementation>>;
    using traits = NeighbourhoodManager_traits<AdaptorMaxOrder>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template<size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;
    //! template<size_t Order, size_t Layer>
    //! using ClusterRefBase_t =
    //!   typename ManagerImplementation::template ClusterRefBase<Order, Layer>;
    //using PairRef_t = ClusterRef_t<2, traits::MaxOrder>;

    // TODO if MaxOrder can be == 1 -> neighbourlist need to be built.
    static_assert(traits::MaxOrder > 1,
                  "ManagerImplementation needs to have an atom list.");

    //! Default constructor
    AdaptorMaxOrder() = delete;

    /**
     * Construct a strict neighbourhood list from a given manager and cut-off
     * radius
     */
    AdaptorMaxOrder(ManagerImplementation& manager, double cutoff);

    //! Copy constructor
    AdaptorMaxOrder(const AdaptorMaxOrder &other) = delete;

    //! Move constructor
    AdaptorMaxOrder(AdaptorMaxOrder &&other) = default;

    //! Destructor
    virtual ~AdaptorMaxOrder() = default;

    //! Copy assignment operator
    AdaptorMaxOrder& operator=(const AdaptorMaxOrder &other) = delete;

    //! Move assignment operator
    AdaptorMaxOrder& operator=(AdaptorMaxOrder &&other) = default;

    //! Update just the adaptor assuming the underlying manager was
    //! updated. this function invokes building either the neighbour list or to
    //! make triplets, quadruplets, etc. depending on the MaxOrder
    void update();

    //! Update the underlying manager as well as the adaptor
    template<class ... Args>
    void update(Args&&... arguments);

    //! Return cutoff radius of the neighbourhood manager
    inline double get_cutoff() const {return this->cutoff;}

    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
				  & counters) const;

    //! Get number of MaxOrder+1 tuplets
    inline size_t get_nb_clusters(size_t cluster_size) const {
      switch (cluster_size) {
      case traits::MaxOrder: {
        return this->neighbours.size();
        break;
      }
      default:
        return this->manager.get_nb_clusters(cluster_size);
        break;
      }
    }

    //! Return number of atoms
    inline size_t get_size() const {
      //! return this->get_nb_clusters(1); TODO: this has to return the original
      //! number of atoms from the underlying manager return
      //! this->neighbours.size();
      return this->manager.get_size();
    }

    //! Return position of atom at given index
    //! (useful for developers)
    inline Vector_ref get_position(const size_t & atom_index) {
      return this->manager.get_position(atom_index);
    }

    //! Return position associated with the given atom object (useful for users)
    inline Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager.get_position(atom.get_index());
    }

    template<size_t Order, size_t Layer>
    inline Vector_ref get_neighbour_position(const ClusterRefBase<Order, Layer>&
                                             /*cluster*/) {
      static_assert(Order > 1,
                    "Only possible for Order > 1.");
      static_assert(Order <= traits::MaxOrder,
                    "this implementation should only work up to MaxOrder.");
      //! Argument is now the same, but implementation
      throw std::runtime_error("should be adapted to Félix's "
                               "new interface using the ClusterRef");
    }

    //! Return the id of an atom at the given index in the neighborhood
    //! of the cluster
    inline int get_cluster_neighbour(const Parent& /*parent*/,
				     size_t index) const {
      return this->manager.get_cluster_neighbour(this->manager, index);
    }

    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefBase<Order, Layer>
				     & cluster,
				     size_t index) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles upto traits::MaxOrder");
      if (Order < traits::MaxOrder-1) {
	return this->manager.get_cluster_neighbour(cluster, index);
      } else {
	auto && offset = this->offsets[cluster.get_cluster_index(Layer)];
	return this->neighbours[offset + index];
      }
    }


    //! Return atom type
    inline int & get_atom_type(const AtomRef_t& atom) {
      return this->manager.get_atom_type(atom.get_index());
    }

    //! Return atom type
    inline const int & get_atom_type(const AtomRef_t& atom) const {
      return this->manager.get_atom_type(atom.get_index());
    }

    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefBase<Order, Layer>
                                   & cluster) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation handles only the respective MaxOrder");
      if (Order < traits::MaxOrder-1) {
	return this->manager.get_cluster_size(cluster);
      } else {
        auto access_index = cluster.get_cluster_index(Layer);
	return nb_neigh[access_index];
      }
    }

    //! Returns the number of neighbors of an atom with the given index
    inline size_t get_cluster_size(const int & atom_index) const {
      return this->manager.get_cluster_size(atom_index);
    }

  protected:
    /**
     * Main function during construction of a neighbourlist.
     *
     * @param atom The atom to add to the list. Because the MaxOrder is
     * increased by one in this adaptor, the Order=MaxOrder
     */
    inline void add_atom(const int atom_index) {
      //! add new atom at this Order
      this->atom_indices.push_back(atom_index);
      //! count that this atom is a new neighbour
      this->nb_neigh.back()++;
      this->offsets.back()++;

      //! make sure that this atom starts with zero lower-Order neighbours
      this->nb_neigh.push_back(0);
      //! update the offsets
      this->offsets.push_back(this->offsets.back() +
                              this->nb_neigh.back());
    }

    //! Create new entry in nb_neigh vector
    inline void add_entry_number_of_neighbours() {
      this->nb_neigh.push_back(0);
    }

    //! Add given atom index as new cluster neighbour
    inline void add_neighbour_of_cluster(const int atom_index) {
      //! add `atom_index` to neighbours
      this->neighbours.push_back(atom_index);
      this->nb_neigh.back()++;
    }

    //! Set the correct offsets for accessing neighbors
    inline void set_offsets() {
      auto n_tuples{nb_neigh.size()};
      this->offsets.reserve(n_tuples);
      this->offsets.resize(1);
      // this->offsets.push_back(0);
      for (size_t i{0}; i<n_tuples; ++i) {
        this->offsets.emplace_back(this->offsets[i] + this->nb_neigh[i]);
      }
    }

    template <size_t Order>
    inline void add_atom(const typename ManagerImplementation::template
                         ClusterRef<Order> & cluster) {
      static_assert(Order <= traits::MaxOrder,
                   "Order too high, not possible to add atom");
      return this->add_atom(cluster.back());
    }

    //! Make a half neighbour list, by construction only Order=1 is supplied.
    void make_half_neighbour_list();

    //! Make a full neighbour list
    void make_full_neighbour_list();

    //! Increase whatever body order is present
    void increase_maxorder();

    ManagerImplementation & manager;

    //! Cutoff radius of manager
    const double cutoff;

    template<size_t Order, bool IsDummy> struct AddOrderLoop;

    //! not necessary any more
    //! // stores AtomRefs to of neighbours for traits::MaxOrder-1-*plets
    // std::vector<AtomRef_t> atom_refs{};

    //! Stores atom indices of current Order
    std::vector<size_t> atom_indices{}; //akin to ilist[]

    //! Stores the number of neighbours for every traits::MaxOrder-1-*plets
    std::vector<size_t> nb_neigh{};

    //! Stores all neighbours of traits::MaxOrder-1-*plets
    std::vector<size_t> neighbours{};

    //! Stores the offsets of traits::MaxOrder-1-*plets for accessing
    //! `neighbours`, from where nb_neigh can be counted
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
  // cutoff goes with the Order.
  template <class ManagerImplementation>
  AdaptorMaxOrder<ManagerImplementation>::
  AdaptorMaxOrder(ManagerImplementation & manager, double cutoff):
    manager{manager},
    cutoff{cutoff},
    atom_indices{},
    nb_neigh{},
    offsets{}

  {
    if (traits::MaxOrder < 1) {
      throw std::runtime_error("No atoms in manager.");
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class ... Args>
  void AdaptorMaxOrder<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <size_t Order, bool IsDummy>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};
    using ClusterRef_t = typename ManagerImplementation::template
      ClusterRef<Order>;

    using NextOrderLoop = AddOrderLoop<Order+1,
				       (Order+1 == OldMaxOrder)>;

    static void loop(ClusterRef_t & cluster,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {

      //! do nothing, if MaxOrder is not reached, except call the next order 
      for (auto next_cluster : cluster) {

        auto & next_cluster_indices
        {std::get<Order>(manager.cluster_indices_container)};
        next_cluster_indices.push_back(next_cluster.get_cluster_indices());

	NextOrderLoop::loop(next_cluster, manager);
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  //! At desired MaxOrder (plus one), here is where the magic happens and the
  //! neighbours of the same order are added as the Order+1.  add check for non
  //! half neighbour list
  template <class ManagerImplementation>
  template <size_t Order>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop<Order, true> {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};

    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;

    // using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    // using IteratorOne_t = typename Manager_t::template iterator<1>;

    using traits = typename AdaptorMaxOrder<ManagerImplementation>::traits;

    static void loop(ClusterRef_t & cluster,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {

      //! get all i_atoms to find neighbours to extend the cluster to the next
      //! order
      auto i_atoms = cluster.get_atom_indices();

      //! vector of existing i_atoms in `cluster` to avoid doubling of atoms in
      //! final list
      std::vector<size_t> current_i_atoms;
      //! a set of new neighbours for the cluster, which will be added to extend
      //! the cluster
      std::set<size_t> current_j_atoms;

      //! access to underlying manager for access to atom pairs
      auto & manager_tmp{cluster.get_manager()};

      // std::cout << " === neigh back "
      //           << i_atoms[0] << " " << i_atoms[1]
      //           << std::endl;

      for (auto atom_index : i_atoms) {
        current_i_atoms.push_back(atom_index);
        size_t access_index = manager.get_cluster_neighbour(manager,
                                                            atom_index);

        //! build a shifted iterator to constuct a ClusterRef<1>
        auto iterator_at_position{manager_tmp.get_iterator_at(access_index)};

        //! ClusterRef<1> as dereference from iterator to get pairs of the
        //! i_atoms
        auto && j_cluster{*iterator_at_position};

        //! collect all possible neighbours of the cluster: collection of all
        //! neighbours of current_i_atoms
        for (auto pair : j_cluster) {
          current_j_atoms.insert(pair.back());
        }
      }

      //! delete existing cluster atoms from list to build additional neighbours
      std::vector<size_t> atoms_to_add{};
      std::set_difference(current_j_atoms.begin(), current_j_atoms.end(),
                          current_i_atoms.begin(), current_i_atoms.end(),
                          std::inserter(atoms_to_add, atoms_to_add.begin()));

      // std::cout << "Neighbours to add to cluster ";
      // for (auto j : atoms_to_add) {
      //   std::cout << j << " ";
      // }
      // std::cout << "number of neighbours to add to add " << atoms_to_add.size();
      // std::cout << std::endl;

      manager.add_entry_number_of_neighbours();
      if (atoms_to_add.size() > 0) {
        for (auto j: atoms_to_add) {
          manager.add_neighbour_of_cluster(j);
        }
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::make_half_neighbour_list() {
    //! Make a half neighbour list (not quite Verlet, because of the missing
    //! skin) according to Tadmor and Miller 'Modeling Materials', algorithm 6.7,
    //! p 323, (needs modification for periodicity).

    //! It results in a <code>strict</code> neighbourlist. It is only a
    //! half-neighbour list.  The most obvious difference is that no 'skin' is
    //! used in conjunction with the cutoff.  This is only necessary, if the
    //! ManagerImplementation, with which this adaptor is initialized does not
    //! have at least already atomic pairs. Inluding ghost neighbours?

    // TODO: add functionality for the shift vector??!

    //! The zeroth order does not have neighbours
    this->nb_neigh.push_back(0);

    unsigned int nneigh_off{0};

    for (auto it=this->manager.begin(); it!=--this->manager.end(); ++it){
      //! Add atom at this order this is just the standard list.
      auto atom_i = *it;
      this->add_atom(atom_i);

      auto jt = it;
      ++jt; //! go to next atom in manager (no self-neighbour)
      for (;jt!=manager.end(); ++jt){
      	auto atom_j = *jt;
      	double distance{(atom_i.get_position() -
			 atom_j.get_position()).norm()};
      	if (distance <= this->cutoff) {
      	  //! Store atom_j in neighbourlist of atom_i
      	  this->nb_neigh.back()++;
      	  nneigh_off += 1;
      	}
      }
      this->offsets.push_back(nneigh_off);
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::make_full_neighbour_list() {
    //! Make a full neighbourlist, whithout fancy linked list or cell. Also
    //! missing are periodic boundary conditions.

    //! The zeroth order does not have neighbours
    this->nb_neigh.push_back(0);

    unsigned int nneigh_off{0};

    for (auto atom_i : this->manager){
      this->add_atom(atom_i);
      for (auto atom_j : this->manager){
      	if(&atom_i != &atom_j) { //! avoid self-neighbouring
	  // TODO: this .get_position() function will always give the
	  // real position, never the shifted one.
      	  double distance{(atom_i.get_position() -
			   atom_j.get_position()).norm()};
      	  if (distance <= this->cutoff) {
      	    //! Store atom_j in neighbourlist of atom_i
	    //! this->atom_refs[1].push_back(atom_j.back());
	    this->add_atom_order_up(atom_j);
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
  void AdaptorMaxOrder<ManagerImplementation>::increase_maxorder() {
    //! Depending on the existing Order, this function increases the
    //! MaxOrder by one by adding all neighbours of the last atom at
    //! the end of the chain as new end of a tuple.
    //! This results in each triplet is only existing once.

    //! Attention, <code>traits::MaxOrder</code> is already increased
    //! in the traits upon construction, therefore the MaxOrder needs
    //! to be larger than 2 (i.e. a NeighbourhoodManager with a
    //! pairlist is present to call this function here.)
    static_assert(traits::MaxOrder > 2, "No neighbourlist present.");

    for (auto atom : this->manager) {
      //! Order 1, Order variable is at 0, atoms, index 0
      using AddOrderLoop = AddOrderLoop<atom.order(),
                                        atom.order() == traits::MaxOrder-1>;
      auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
      atom_cluster_indices.push_back(atom.get_cluster_indices());
      AddOrderLoop::loop(atom, *this);
    }

    //! correct the offsets for the new cluster order
    this->set_offsets();
    //! add correct cluster_indices for the highest order
    auto & max_cluster_indices
    {std::get<traits::MaxOrder-1>(this->cluster_indices_container)};
    max_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::update() {
    //! initialise the list if it does not exist
    // TODO
    // initialise the neighbourlist

    //! -1 because the traits' MaxOrder is already increased
    if (traits::MaxOrder-1 == 1) {
      //! Make half neighbour list (strict?)
      //! initialise the neighbourlist
      for (int i{0}; i < traits::MaxOrder; ++i) {
	//! this->atom_refs.clear();
	this->nb_neigh.resize(0);
	this->offsets.resize(0);
      }
      this->make_half_neighbour_list();
      // this->make_full_neighbour_list(); // no frills, full neighbourlist
    } else {
      //! Make triplets/quadruplets/etc. based on existing
      //! neighbourlist
      //! Templated function?
      this->increase_maxorder();
    }
  }

  /* ---------------------------------------------------------------------- */
  template<class ManagerImplementation>
  template<size_t Order>
  inline size_t AdaptorMaxOrder<ManagerImplementation>::
  get_offset_impl(const std::array<size_t, Order> & counters) const {

    static_assert(Order < traits::MaxOrder,
                  "this implementation handles only the respective MaxOrder");
    /**
     * Order accessor: 0 - atoms
     *                 1 - pairs
     *                 2 - triplets
     *                 etc.
     * Order is determined by the ClusterRef building iterator, not by the Order
     * of the built iterator
     */

    if (Order == traits::MaxOrder-1) {
      /**
       * Counters as an array to call parent offset multiplet. This can then be
       * used to access the actual offset for the next Order here.
       */

      auto i{this->manager.get_offset_impl(counters)};
      auto j{counters[Order-1]};
      auto tuple_index{i+j};
      auto main_offset{this->offsets[tuple_index]};
      return main_offset;
    } else {
      /**
       * If not accessible at this level, call lower Order offsets from lower
       * order manager(s).
       */
      return this->manager.get_offset_impl(counters);
    }
  }
}  // rascal

#endif /* ADAPTOR_MAXLEVEL_H */

// TODO: The construction of triplets is fine, but they occur multiple times. We
// probably need to check for increasing atomic index to get rid of
// duplicates. But this is in general a design decision, if we want full/half neighbour list and full/half/whatever triplets and quadruplets
