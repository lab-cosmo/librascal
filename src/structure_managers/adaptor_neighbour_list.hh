/**
 * file   adaptor_neighbour_list.hh
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



#ifndef ADAPTOR_NEIGHBOUR_LIST_H
#define ADAPTOR_NEIGHBOUR_LIST_H

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"
#include "lattice.hh"
#include "basic_types.hh"

#include <typeinfo>
#include <set>
#include <vector>

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
     * Mesh bounding coordinates iterator for easy access to the cornes for
     * evaluating the multipliers necessary to build periodic images.
     */
    template <size_t Dim>
    class MeshBounds {
    public:
      //! constructor
      MeshBounds(const std::array<double, 2*Dim> & extent)
        : extent{extent} {};
      //! copy constructor
      MeshBounds(const MeshBounds & other) = default;
      //! assignment operator
      MeshBounds & operator=(const MeshBounds & other) = default;
      ~MeshBounds() = default;

      //! iterators over `` dereferences to mesh bound coordinate
      class iterator
      {
      public:
        using value_type = std::array<double, Dim>; //!< stl conformance
        using const_value_type = const value_type; //!< stl conformance
        using pointer = value_type*; //!< stl conformance
        using iterator_category = std::forward_iterator_tag;//!<stl conformance

        //! constructor
        iterator(const MeshBounds & mesh_bounds, bool begin=true)
          : mesh_bounds{mesh_bounds},
            index{begin? 0: mesh_bounds.size()} {}

        ~iterator() {};
        //! dereferencing
        value_type operator*() const {
          std::array<double, Dim> retval{{0}};
          constexpr int size{2};
          for (size_t i{0}; i < Dim; ++i) {
            int idx = (this->index/ipow(size,i))%size * Dim + i;
            retval[i] = this->mesh_bounds.extent[idx];
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
        const MeshBounds & mesh_bounds;
        //!< index of currect pointed-to pixel
        size_t index;
      };
      //! stl conformance
      inline iterator begin() const {return iterator(*this);}
      //! stl conformance
      inline iterator end() const {return iterator(*this, false);}
      //! stl conformance
      inline size_t size() const {return ipow(2, Dim);}
    protected:
      const std::array<double, 2*Dim> extent; //!< repetitions in each dimension
    };

    /* ---------------------------------------------------------------------- */
    //! get dimension dependent neighbour indices
    template<size_t Dim, class Container_t>
    std::vector<size_t> get_neighbours(const int current_atom_index,
                                       const std::array<int, Dim> & ccoord,
                                       const Container_t & boxes) {
      std::vector<size_t> neighbours;
      for (auto && s: Stencil<Dim>{ccoord}) {
        for (const auto & neigh : boxes[s]) {
          //! avoid adding the current i atom to the neighbour list
          if (neigh != current_atom_index) {
            neighbours.push_back(neigh);
          }
        }
      }
      return neighbours;
    }

    /* ---------------------------------------------------------------------- */
    //! get the cell index for a position
    template<class Vector_t>
    decltype(auto) get_box_index(const Vector_t & position,
                                 const double & rc) {
      // const std::array<int, Dim> nmax) {

      auto constexpr dimension{Vector_t::SizeAtCompileTime};
      // static_assert(dimension == Dim,
      //               "Discrepancy between position and boxgrid dimension");

      std::array<int, dimension> nidx{};
      for(auto dim{0}; dim < dimension; ++dim) {
        auto val = position(dim);
        nidx[dim] = int(std::floor(val / rc));
        // std::cout << "dim, pos " << dim << " " << val;
        // std::cout << " val/rc " << val/rc << " " << int(std::floor(val / rc)) << std::endl;
        // nidx[dim] = std::min(nidx[dim], nmax[dim]-1);
        // nidx[dim] = std::max(nidx[dim], 0);
      }
      return nidx;
    }

    /* ---------------------------------------------------------------------- */
    //! get the linear index of a voxel in a given grid
    template <size_t Dim>
    constexpr Dim_t get_index(const std::array<int, Dim> & sizes,
                              const std::array<int, Dim> & ccoord) {
      Dim_t retval{0};
      Dim_t factor{1};
      for (Dim_t i = Dim-1; i >=0; --i) {
        retval += ccoord[i]*factor;
        if (i != 0) {
          factor *= sizes[i];
        }
      }
      return retval;
    }

    /* ---------------------------------------------------------------------- */
    //! get the dim-index array from a linear index
    template <size_t Dim>
    constexpr std::array<int, Dim>
    get_ccoord(const std::array<int, Dim> & sizes,
               const std::array<int, Dim> & origin, int index) {
      std::array<int, Dim> retval{{0}};
      int factor{1};
      for (size_t i = Dim-1; i >=0; --i) {
        retval[i] = index/factor%sizes[i] + origin[i];
        if (i != 0 ) {
          factor *= sizes[i];
        }
      }
      return retval;
    }

    /* ---------------------------------------------------------------------- */
    //! test if position inside
    template <int Dim>
    bool position_in_bounds(const Eigen::Matrix<double, Dim, 1> & min,
                            const Eigen::Matrix<double, Dim, 1> & max,
                            const Eigen::Matrix<double, Dim, 1> & pos) {

      auto pos_lower  = pos.array() - min.array();
      auto pos_greater  = pos.array() - max.array();

      //! check if shifted position inside mesh
      auto f_lt = (pos_lower.array() > 0.).all();
      auto f_gt = (pos_greater.array() < 0.).all();

      if (f_lt and f_gt) {
        return true;
      } else {
        return false;
      }
    }

    /* ---------------------------------------------------------------------- */
    /**
     * storage for cell coordinates of atoms depending on the number of
     * dimensions
     */
    template<int Dim>
    class IndexContainer
    {
    public:
      //! Default constructor
      IndexContainer() = delete;

      //! Constructor with size
      IndexContainer(const std::array<int, Dim> & nboxes)
        : nboxes{nboxes} {
        auto ntot = std::accumulate(nboxes.begin(), nboxes.end(),
                                    1, std::multiplies<int>());
        data.resize(ntot);
      }

      //! Copy constructor
      IndexContainer(const IndexContainer &other) = delete;

      //! Move constructor
      IndexContainer(IndexContainer &&other) = delete;

      //! Destructor
      ~IndexContainer(){};

      //! Copy assignment operator
      IndexContainer& operator=(const IndexContainer &other) = delete;

      //! Move assignment operator
      IndexContainer& operator=(IndexContainer &&other) = default;

      std::vector<int> & operator[](const std::array<int, Dim>& ccoord) {
        auto index = get_index(this->nboxes, ccoord);
        return data[index];
      }

      const std::vector<int> & operator[](const std::array<int,
                                          Dim>& ccoord) const {
        auto index = get_index(this->nboxes, ccoord);
        return this->data[index];
      }

    protected:
      //! a vector of atom indices for every box
      std::vector<std::vector<int>> data{};
      //! number of boxes in each dimension
      std::array<int, Dim> nboxes{};
    private:
    };

  }  // internal


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
    using AtomTypes_ref = Eigen::Map<Eigen::Matrix<int, 1, Eigen::Dynamic>>;

    static_assert(traits::MaxOrder == 2,
                  "ManagerImplementation needs an atom list "
                  " and can only build a neighbour list (pairs).");

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

    //! ghost types are only available for MaxOrder=2
    inline int & get_ghost_type(const size_t & atom_index) {
      auto p = this->get_ghost_types();
      return p(atom_index);
    }

    inline AtomTypes_ref get_ghost_types() {
      AtomTypes_ref val(this->ghost_types.data(), 1, this->ghost_types.size());
      return val;
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

      return this->get_position(cluster.back());
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

    //! Returns atom type given an atom index, also works for ghost atoms
    inline int get_atom_type(const size_t & atom_index) {
      if (atom_index < this->n_i_atoms) {
        return this->manager.get_atom_type(atom_index);
      } else {
        return this->get_ghost_type(atom_index - this->n_i_atoms);
      }
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
      //! if current adaptor has built pair list
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
    inline void add_ghost_atom(const int & atom_index,
                               const Vector_t & position,
                               const int & atom_type) {
      this->atom_indices.push_back(atom_index);
      this->ghost_types.push_back(atom_type);
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

    //! ghost atom positions
    std::vector<double> ghost_positions{};

    //! ghost atom type
    std::vector<int> ghost_types{};
  private:
  };

  /* ---------------------------------------------------------------------- */
  //! Construction of the pair list manager

  template <class ManagerImplementation>
  AdaptorNeighbourList<ManagerImplementation>::
  AdaptorNeighbourList(ManagerImplementation & manager, double cutoff):
    manager{manager},
    cutoff{cutoff},
    atom_indices{},
    nb_neigh{},
    offsets{}
  {
    if (traits::MaxOrder < 1) {
      throw std::runtime_error("No atom list in manager.");
    }
    //! reset number of center atoms and initialize number of ghosts
    n_i_atoms = this->manager.get_size();
    n_j_atoms = 0;
  }


  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorNeighbourList<ManagerImplementation>::update() {
    /**
     * build a neighbour list based on atomic positions, types and indices, in
     * the following the needed data structures are initialized
     */
    this->nb_neigh.resize(0);
    this->offsets.resize(0);
    this->neighbours.resize(0);
    this->ghost_types.resize(0);
    this->make_full_neighbour_list();
  }
  /* ---------------------------------------------------------------------- */
  /**
   * Builds full neighbour list. Triclinicity is accounted for. The general idea
   * is to anchor a mesh at the origin of the supplied cell (assuming it is at
   * the origin). Then the mesh is extended into space until it is as big as the
   * maximum cell coordinate plus one cutoff in each direction. This mesh has
   * boxes of size ``cutoff``. Depending on the periodicity of the mesh, ghost
   * atoms are added by shifting all i-atoms by the cell vectors corresponding
   * to the desired periodicity. All i-atoms and the ghost atoms are then sorted
   * into the respective boxes of the cartesian mesh and a stencil anchored at
   * the box of the i-atoms is used to build the atom neighbourhoods from the 9
   * (2d) or 27 (3d) boxes, which is checked for the neighbourhood. Correct
   * periodicity is ensured by the placement of the ghost atoms. The resulting
   * neighbourlist is full and not strict.
   */
  template <class ManagerImplementation>
  void AdaptorNeighbourList<ManagerImplementation>::make_full_neighbour_list() {
    using Vector_t = Eigen::Matrix<double, traits::Dim, 1>;

    //! short hands for variable
    const auto dim{traits::Dim};

    auto cell{this->manager.get_cell()};
    double cutoff{this->cutoff};
    std::array<int, dim> nboxes_per_dim;

    //! vector for storing the atom indices of each box
    std::vector<std::vector<int>> atoms_in_box{};

    /**
     * minimum/maximum coordinate of mesh for neighbour list; depends on cell
     * triclinicity and cutoff, coordinates of the mesh are relative to the
     * origin of the given cell.
     */
    Vector_t mesh_min(dim);
    Vector_t mesh_max(dim);

    mesh_min.setZero();
    mesh_max.setZero();

    /**
     * max and min multipliers for number of cells in mesh per dimension in
     * units of cell vectors
     */
    std::array<int, dim> m_min;
    std::array<int, dim> m_max;

    /**
     * Mesh related stuff for neighbour boxes. Calculate min and max of the mesh
     * in cartesian coordinates and relative to the cell origin.  mesh_min is
     * the origin of the mesh; mesh_max is the maximum coordinate of the mesh;
     * nboxes_per_dim is the number of mesh boxes in each dimension, not to be
     * confused with the number of cells to ensure periodicity
     */
    for (auto i{0}; i < dim; ++i) {
      auto min_coord = std::min(0., cell.row(i).minCoeff());
      auto max_coord = std::max(0., cell.row(i).maxCoeff());

      /**
       * minimum is given by -cutoff and a delta to avoid ambiguity during cell
       * sorting of atom position e.g. at x = (0,0,0).
       */
      auto epsilon = 0.25*cutoff;
      mesh_min[i] = min_coord - cutoff - epsilon;
      auto lmesh = std::fabs(mesh_min[i]) + max_coord + cutoff;
      int n = std::ceil(lmesh / cutoff);
      auto lmax = n * cutoff - std::fabs(mesh_min[i]);
      mesh_max[i] = lmax;
      nboxes_per_dim[i] = n;
    }

    /**
     * Periodicity related multipliers. Now the mesh coordinates are calculated
     * in units of cell vectors. m_min and m_max give the number of repetitions
     * of the cell in each cell vector direction
     */
    int ncorners = internal::ipow(2, dim);
    Eigen::MatrixXd xpos(dim, ncorners);
    std::array<double, dim*2> mesh_bounds{};
    for (auto i{0}; i<dim; ++i) {
      mesh_bounds[i] = mesh_min[i];
      mesh_bounds[i+dim] = mesh_max[i];
    }
    // TODO: find way to initialize eigen column with std::array
    int n{0};
    for (auto && coord : internal::MeshBounds<dim>{mesh_bounds}) {
      for (auto i{0}; i<dim; ++i) {
        xpos(i,n) = coord[i];
      }
      n++;
    }
    //! solve for all multipliers
    auto multiplicator = cell.ldlt().solve(xpos);
    auto xmin = multiplicator.rowwise().minCoeff();
    auto xmax = multiplicator.rowwise().maxCoeff();

    for (auto i{0}; i < dim; ++i) {
      //! +/- 1 because of the "zero" cell, the cell itself
      m_min[i] = std::floor(xmin(i)) - 1;
      m_max[i] = std::ceil(xmax(i)) + 1;
    }

    /**
     * TODO possible future optimization for cells large triclinicity: use
     * triclinic coordinates and explicitly check for the 'skin' around the cell
     * and rotate the cell to have the lower triangular form
     */
    std::array<int, dim> periodic_max{};
    std::array<int, dim> periodic_min{};
    std::array<int, dim> repetitions{};
    auto periodicity = this->manager.get_periodic_boundary_conditions();
    size_t ntot{1};

    //! calculate number of actual repetitions of cell, depending on periodicity
    for (auto i{0}; i < dim; ++i) {
      if (periodicity[i]) {
        periodic_max[i] = m_max[i];
        periodic_min[i] = m_min[i];
      } else {
        periodic_max[i] = 0;
        periodic_min[i] = 0;
      }
      auto nrep_in_dim = -periodic_min[i] + periodic_max[i] + 1;
      repetitions[i] = nrep_in_dim;
      ntot *= nrep_in_dim;
    }

    //! generate ghost atom indices and position
    for (auto atom : this->get_manager()) {
      auto pos = atom.get_position();
      auto atom_type = atom.get_atom_type();

      for (auto && p_image : internal::PeriodicImages<dim>
        {periodic_min, repetitions, ntot}) {

        int ncheck{0};
        for (auto i{0}; i<dim; ++i) ncheck += std::abs(p_image[i]);

        //! exclude cell itself
        if(ncheck > 0) {
          Vector_t pos_ghost = pos;

          for (auto i{0}; i< dim; ++i) {
            pos_ghost += cell.col(i) * p_image[i];
          }

          auto flag_inside = internal::position_in_bounds(mesh_min,
                                                          mesh_max,
                                                          pos_ghost);
          if (flag_inside) {
            //! next atom index is size, since start is at index = 0
            auto new_atom_index = this->get_size_with_ghosts();
            this->add_ghost_atom(new_atom_index, pos_ghost, atom_type);
          }
        }
      }
    }

    //! neighbour boxes
    internal::IndexContainer<dim> atom_id_cell{nboxes_per_dim};

    // i-atoms sorting into boxes
    for (size_t i{0}; i < this->n_i_atoms; ++i) {
      Vector_t pos = this->get_position(i);
      Vector_t dpos = pos - mesh_min;
      auto idx = internal::get_box_index(dpos, cutoff);
      atom_id_cell[idx].push_back(i);
    }

    //! ghost atoms sorting into boxes
    for (size_t i{0}; i < this->n_j_atoms; ++i) {
      Vector_t ghost_pos = this->get_ghost_position(i);
      Vector_t dpos = ghost_pos - mesh_min;
      auto idx  = internal::get_box_index(dpos, cutoff);
      auto ghost_atom_index = i + this->n_i_atoms;
      atom_id_cell[idx].push_back(ghost_atom_index);
    }

    //! go through atoms and build neighbour list
    int offset{0};
    for (size_t i{0}; i < this->n_i_atoms; ++i) {
      int nneigh{0};
      Vector_t pos = this->get_position(i);
      Vector_t dpos = pos - mesh_min;
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
