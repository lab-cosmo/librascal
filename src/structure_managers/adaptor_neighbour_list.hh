/**
 * file   adaptor_increase_maxorder.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   04 Oct 2018
 *
 * @brief implements an adaptor for structure_managers, which
 * creates a full or half neighbourlist if there is none
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

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"
#include "lattice.hh"
#include "basic_types.hh"

#include <typeinfo>
#include <set>
#include <vector>

#ifndef ADAPTOR_NEIGHBOUR_LIST_H
#define ADAPTOR_NEIGHBOUR_LIST_H

namespace rascal {
  /**
   * Forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorNeighbourList;

  /**
   * Specialisation of traits for increase <code>MaxOrder</code> adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorNeighbourList<ManagerImplementation>> {

    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
      ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    //! New MaxOrder upon construction!
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder+1};
    //! New Layer
    using LayerByOrder =
      typename LayerExtender<MaxOrder,
                             typename
                             ManagerImplementation::traits::LayerByOrder>::type;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Adaptor that increases the MaxOrder of an existing StructureManager. This
   * means, if the manager does not have a neighbourlist, it is created, if it
   * exists, triplets, quadruplets, etc. lists are created.
   */
  template <class ManagerImplementation>
  class AdaptorNeighbourList: public
  StructureManager<AdaptorNeighbourList<ManagerImplementation>>
  {
  public:
    using Base = StructureManager<AdaptorNeighbourList<ManagerImplementation>>;

    using Parent =
      StructureManager<AdaptorNeighbourList<ManagerImplementation>>;
    using traits = StructureManager_traits<AdaptorNeighbourList>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template<size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using Positions_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
                                                   Eigen::Dynamic>>;

    static_assert(traits::MaxOrder == 1,
                  "ManagerImplementation needs an atom list "
                  " and can only build a neighbour list.");

    //! Default constructor
    AdaptorNeighbourList() = delete;

    /**
     * Constructs a full neighbourhood list from a given manager and cut-off
     * radius or extends an existing neighbourlist to the next order
     */
    AdaptorNeighbourList(ManagerImplementation & manager, double cutoff);

    //! Copy constructor
    AdaptorNeighbourList(const AdaptorNeighbourList & other) = delete;

    //! Move constructor
    AdaptorNeighbourList(AdaptorNeighbourList && other) = default;

    //! Destructor
    virtual ~AdaptorNeighbourList() = default;

    //! Copy assignment operator
    AdaptorNeighbourList & operator=(const AdaptorNeighbourList & other) = delete;

    //! Move assignment operator
    AdaptorNeighbourList & operator=(AdaptorNeighbourList && other) = default;

    /**
     * Updates just the adaptor assuming the underlying manager was
     * updated. this function invokes building either the neighbour list or to
     * make triplets, quadruplets, etc. depending on the MaxOrder
     */
    void update();

    //! Updates the underlying manager as well as the adaptor
    template<class ... Args>
    void update(Args&&... arguments);

    //! Returns cutoff radius of the neighbourhood manager
    inline double get_cutoff() const {return this->cutoff;}

    /**
     * Returns the linear indices of the clusters (whose atom indices are stored
     * in counters). For example when counters is just the list of atoms, it
     * returns the index of each atom. If counters is a list of pairs of indices
     * (i.e. specifying pairs), for each pair of indices i,j it returns the
     * number entries in the list of pairs before i,j appears.
     */
    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
				  & counters) const;

    //! Returns the number of clusters of size cluster_size
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

    //! Returns number of clusters of the original manager
    inline size_t get_size() const {
      return this->manager.get_size();
    }

    //! total number of atoms used for neighbour list, including ghosts
    inline size_t get_size_with_ghosts() const{
      return this->n_i_atoms+this->n_j_atoms;
    }

    //! Returns position of an atom with index atom_index
    inline Vector_ref get_position(const size_t & atom_index) {
      if (atom_index < n_i_atoms) {
        return this->manager.get_position(atom_index);
      } else {
        return this->get_ghost_position(atom_index - this->n_i_atoms);
      }
    }

    //! ghost positions are only available for MaxOrder == 2
    inline Vector_ref get_ghost_position(const size_t & atom_index) {
      auto p = this->get_ghost_positions();
      auto * xval{p.col(atom_index).data()};
      return Vector_ref(xval);
    }

    inline Positions_ref get_ghost_positions() {
      return Positions_ref(this->ghost_positions.data(), traits::Dim,
                           this->ghost_positions.size()/traits::Dim);
    }

    //! Returns position of the given atom object (useful for users)
    inline Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager.get_position(atom.get_index());
    }

    template<size_t Order, size_t Layer>
    inline Vector_ref get_neighbour_position(const ClusterRefKey<Order, Layer>
                                             & cluster) {
      static_assert(Order > 1,
                    "Only possible for Order > 1.");
      static_assert(Order <= traits::MaxOrder,
                    "this implementation should only work up to MaxOrder.");

      // if (Order == 2) {
      return this->get_position(cluster.back());
      // } else {
      //return this->get_position(cluster.back());//manager.get_neighbour_position(cluster);
      //}
    }

    /**
     * Returns the id of the index-th (neighbour) atom of the cluster that is
     * the full structure/atoms object, i.e. simply the id of the index-th atom
     */
    inline int get_cluster_neighbour(const Parent& /*parent*/,
				     size_t index) const {
      return this->manager.get_cluster_neighbour(this->manager, index);
    }

    //! Returns the id of the index-th neighbour atom of a given cluster
    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
				     & cluster,
				     size_t index) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles up to traits::MaxOrder");
      if (Order < traits::MaxOrder-1) {
	return this->manager.get_cluster_neighbour(cluster, index);
      } else {
	auto && offset = this->offsets[cluster.get_cluster_index(Layer)];
	return this->neighbours[offset + index];
      }
    }

    //! Returns atom type given an atom object AtomRef
    inline int & get_atom_type(const AtomRef_t& atom) {
      return this->manager.get_atom_type(atom.get_index());
    }

    //! Returns a constant atom type given an atom object AtomRef
    inline const int & get_atom_type(const AtomRef_t& atom) const {
      return this->manager.get_atom_type(atom.get_index());
    }

    //! Returns the number of neighbors of a given cluster
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
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
      if (traits::MaxOrder == 2) {
        return this->nb_neigh[atom_index];
      } else {
        return this->manager.get_cluster_size(atom_index);
      }
    }

  protected:
    /**
     * Main function during construction of a neighbourlist.
     *
     * @param atom The atom to add to the list. Because the MaxOrder is
     * increased by one in this adaptor, the Order=MaxOrder
     */
    inline void add_atom(const int atom_index) {
      //! adds new atom at this Order
      this->atom_indices.push_back(atom_index);
      //! increases the number of neighbours
      this->nb_neigh.back()++;
      //! increases the offsets
      this->offsets.back()++;

      /**
       * extends the list containing the number of neighbours with a new 0 entry
       * for the added atom
       */
      this->nb_neigh.push_back(0);

      /**
       * extends the list containing the offsets and sets it with the number of
       * neighbours plus the offsets of the last atom
       */
      this->offsets.push_back(this->offsets.back() +
                              this->nb_neigh.back());
    }

    /**
     * This function, including the storage of ghost atom positions is
     * necessary, because the underlying manager is not known at this
     * layer. Therefore we can not add positions to the existing array, but have
     * to add positions to a ghost array. This also means, that the get_position
     * function will need to branch, depending on the atom_index > n_i_atoms and
     * offset with n_j_atoms to access ghost positions.
     */
    inline void add_ghost_atom(const int atom_index, const Vector_t position) {
      this->atom_indices.push_back(atom_index);
      for (auto dim{0}; dim < traits::Dim; ++dim) {
        this->ghost_positions.push_back(position(dim));
      }
      this->n_j_atoms++;
    }

    //! Extends the list containing the number of neighbours with a 0
    inline void add_entry_number_of_neighbours() {
      this->nb_neigh.push_back(0);
    }

    //! Adds a given atom index as new cluster neighbour
    inline void add_neighbour_of_cluster(const int atom_index) {
      //! adds `atom_index` to neighbours
      this->neighbours.push_back(atom_index);
      //! increases the number of neighbours
      this->nb_neigh.back()++;
    }

    //! Sets the correct offsets for accessing neighbours
    inline void set_offsets() {
      auto n_tuples{nb_neigh.size()};
      this->offsets.reserve(n_tuples);
      this->offsets.resize(1);
      // this->offsets.push_back(0);
      for (size_t i{0}; i < n_tuples; ++i) {
        this->offsets.emplace_back(this->offsets[i] + this->nb_neigh[i]);
      }
    }

    /**
     * Interface of the add_atom function that adds the last atom in a given
     * cluster
     */
    template <size_t Order>
    inline void add_atom(const typename ManagerImplementation::template
                         ClusterRef<Order> & cluster) {
      static_assert(Order <= traits::MaxOrder,
                    "Order too high, not possible to add atom");
      return this->add_atom(cluster.back());
    }

    //! full neighbour list with cell algorithm if Order==1
    void make_full_neighbour_list();

    ManagerImplementation & manager;

    //! Cutoff radius of manager
    const double cutoff;

    template<size_t Order, bool IsDummy> struct AddOrderLoop;

    /**
     * Compile time decision, if a new neighbour list is built or if an existing
     * one is extended to the next Order.
     */
    template <size_t Order, bool IsDummy> struct IncreaseMaxOrder;

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

    //! number of i atoms, i.e. centers
    size_t n_i_atoms{};
    /**
     * number of ghost atoms (given by periodicity) filled during full
     * neighbourlist build
     */
    size_t n_j_atoms{};

    //! ghost positions
    std::vector<double> ghost_positions{};
  private:
  };

  namespace internal {

    /* ---------------------------------------------------------------------- */
    template <typename R, typename I>
    constexpr R ipow(R base, I exponent) {
      static_assert(std::is_integral<I>::value, "Type must be integer");
      R retval{1};
      for (I i = 0; i < exponent; ++i) {
        retval *= base;
      }
      return retval;
    }

    /* ---------------------------------------------------------------------- */
    /**
     * stencil iterator for simple, dimension dependent stencils to access the
     * neighbouring boxes of the cell algorithm
     */
    template <size_t Dim>
    class Stencil {
    public:
      //! constructor
      Stencil(const std::array<int, Dim> & origin)
        : origin{origin}{};
      //! copy constructor
      Stencil(const Stencil & other) = default;
      //! assignment operator
      Stencil & operator=(const Stencil & other) = default;
      ~Stencil() = default;

      //! iterators over `` dereferences to cell coordinates
      class iterator
      {
      public:
        using value_type = std::array<int, Dim>; //!< stl conformance
        using const_value_type = const value_type; //!< stl conformance
        using pointer = value_type*; //!< stl conformance
        using iterator_category = std::forward_iterator_tag;//!<stl conformance

        //! constructor
        iterator(const Stencil & stencil, bool begin=true)
          : stencil{stencil}, index{begin? 0: stencil.size()} {}

        ~iterator() {};
        //! dereferencing
        value_type operator*() const {
          constexpr int size{3};
          std::array<int, Dim> retval{{0}};
          int factor{1};
          for (int i = Dim-1; i >=0; --i) {
            //! -1 for offset of stencil
            retval[i] = this->index/factor%size + this->stencil.origin[i] - 1;
            if (i != 0 ) {
              factor *= size;
            }
          }
          return retval;
        };
        //! pre-increment
        iterator & operator++() {this->index++; return *this;}
        //! inequality
        inline bool operator!=(const iterator & other) const {
          return this->index != other.index;
        };

      protected:
        //!< ref to stencils in cell
        const Stencil & stencil;
        //!< index of currect pointed-to pixel
        size_t index;
      };
      //! stl conformance
      inline iterator begin() const {return iterator(*this);}
      //! stl conformance
      inline iterator end() const {return iterator(*this, false);}
      //! stl conformance
      inline size_t size() const {return ipow(3, Dim);}
    protected:
      const std::array<int, Dim> origin; //!< locations of this domain
    };

    /* ---------------------------------------------------------------------- */
    /**
     * Periodic image iterator for easy access to how many images have to be
     * added for ghost atoms.
     */
    template <size_t Dim>
    class PeriodicImages {
    public:
      //! constructor
      PeriodicImages(const std::array<int, Dim> & origin,
                     const std::array<int, Dim> & nrepetitions,
                     const size_t & ntot)
        : origin{origin}, nrepetitions{nrepetitions}, ntot{ntot} {};
      //! copy constructor
      PeriodicImages(const PeriodicImages & other) = default;
      //! assignment operator
      PeriodicImages & operator=(const PeriodicImages & other) = default;
      ~PeriodicImages() = default;

      //! iterators over `` dereferences to cell coordinates
      class iterator
      {
      public:
        using value_type = std::array<int, Dim>; //!< stl conformance
        using const_value_type = const value_type; //!< stl conformance
        using pointer = value_type*; //!< stl conformance
        using iterator_category = std::forward_iterator_tag;//!<stl conformance

        //! constructor
        iterator(const PeriodicImages & periodic_images, bool begin=true)
          : periodic_images{periodic_images},
            index{begin? 0: periodic_images.size()} {}

        ~iterator() {};
        //! dereferencing
        value_type operator*() const {
          std::array<int, Dim> retval{{0}};
          int factor{1};
          for (int i = Dim-1; i >=0; --i) {
            retval[i] = this->index/factor%this->periodic_images.nrepetitions[i]
              + this->periodic_images.origin[i];
            if (i != 0 ) {
              factor *= this->periodic_images.nrepetitions[i];
            }
          }
          return retval;
        };
        //! pre-increment
        iterator & operator++() {this->index++; return *this;}
        //! inequality
        inline bool operator!=(const iterator & other) const {
          return this->index != other.index;
        };

      protected:
        //!< ref to periodic images
        const PeriodicImages & periodic_images;
        //!< index of currect pointed-to pixel
        size_t index;
      };
      //! stl conformance
      inline iterator begin() const {return iterator(*this);}
      //! stl conformance
      inline iterator end() const {return iterator(*this, false);}
      //! stl conformance
      inline size_t size() const {return this->ntot;}
    protected:
      const std::array<int, Dim> origin; //!< minimum repetitions
      const std::array<int, Dim> nrepetitions; //!< repetitions in each dimension
      const size_t ntot;
    };

    /* ---------------------------------------------------------------------- */
    /**
     * Mesh bounding coordinates iterator for easy access to the corne;
      auto idx = internal::get_box_index(dpos, cutoff);
      auto current_j_atoms = internal::get_neighbours(i, idx, atom_id_cell);

      for (auto j : current_j_atoms) {
        this->neighbours.push_back(j);
        nneigh++;
      }
      this->nb_neigh.push_back(nneigh);
      this->offsets.push_back(offset);
      offset += nneigh;
    }

    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    atom_cluster_indices.fill_sequence();
    pair_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorNeighbourList<ManagerImplementation>::update() {
    /**
     * Standard case, increase an existing neighbour list or triplet list to a
     * higher Order
     */
    manager_max->nb_neigh.resize(0);
    manager_max->offsets.resize(0);
    manager_max->neighbours.resize(0);
    manager_max->make_full_neighbour_list();
  }

  /* ---------------------------------------------------------------------- */
  /* Returns the linear indices of the clusters (whose atom indices
   * are stored in counters). For example when counters is just the list
   * of atoms, it returns the index of each atom. If counters is a list of pairs
   * of indices (i.e. specifying pairs), for each pair of indices i,j it returns
   * the number entries in the list of pairs before i,j appears.
   */
  template<class ManagerImplementation>
  template<size_t Order>
  inline size_t AdaptorNeighbourList<ManagerImplementation>::
  get_offset_impl(const std::array<size_t, Order> & counters) const {

    static_assert(Order < traits::MaxOrder,
                  "this implementation handles only up to "
                  "the respective MaxOrder");
    /**
     * Order accessor: 0 - atoms
     *                 1 - pairs
     *                 2 - triplets
     *                 etc.
     * Order is determined by the ClusterRef building iterator, not by the Order
     * of the built iterator
     */
    if (Order == 1) {
      return this->offsets[counters.front()];
    } else if (Order == 2) {
      /**
       * Counters as an array to call parent offset multiplet. This can then be
       * used to access the actual offset for the next Order here.
       */
      auto i{this->manager.get_offset_impl(counters)};
      auto j{counters[Order-1]};
      auto tuple_index{i+j};
      auto main_offset{this->offsets[tuple_index]};
      return main_offset;
    }
  }
}  // rascal

#endif /* ADAPTOR_NEIGHBOUR_LIST_H */
