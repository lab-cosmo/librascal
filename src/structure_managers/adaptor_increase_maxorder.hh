/**
 * file   adaptor_increase_maxorder.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   19 Jun 2018
 *
 * @brief implements an adaptor for structure_managers, which
 * creates a full and half neighbourlist if there is none and
 * triplets/quadruplets, etc. if existent.
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

#ifndef ADAPTOR_MAXORDER_H
#define ADAPTOR_MAXORDER_H

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
  struct StructureManager_traits<AdaptorMaxOrder<ManagerImplementation>> {

    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
      ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    //! New MaxOrder upon construction!
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder+1};
    //! New Layer
    //! TODO: Is this the correct way to initialize the increased order?
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
  class AdaptorMaxOrder: public
  StructureManager<AdaptorMaxOrder<ManagerImplementation>>
  {
  public:
    using Base = StructureManager<AdaptorMaxOrder<ManagerImplementation>>;

    using Parent =
      StructureManager<AdaptorMaxOrder<ManagerImplementation>>;
    using traits = StructureManager_traits<AdaptorMaxOrder>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template<size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using Positions_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
                                                   Eigen::Dynamic>>;

    static_assert(traits::MaxOrder > 1,
                  "ManagerImplementation needs an atom list.");

    //! Default constructor
    AdaptorMaxOrder() = delete;

    /**
     * Constructs a full neighbourhood list from a given manager and cut-off
     * radius or extends an existing neighbourlist to the next order
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
      static_assert(Order < traits::MaxOrder,
                    "this implementation should only work up to MaxOrder.");

      if (Order == 2) {
        return this->get_position(cluster.back());
      } else {
        return this->manager.get_neighbour_position(cluster, index);
      }
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
    //! stencil iterator for simple, dimension dependent stencils
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
        const Stencil& stencil;
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
    //! get the cell index or a position
    template<class Vector_t, size_t Dim>
    decltype(auto) get_box_index(const Vector_t & position,
                                 const double & rc,
                                 const std::array<int, Dim> nmax) {

      auto constexpr dimension{Vector_t::SizeAtCompileTime};
      static_assert(dimension == Dim,
                    "Discrepancy between position and boxgrid dimension");

      std::array<int, dimension> nidx{};
      for(auto dim{0}; dim < dimension; ++dim) {
        nidx[dim] = static_cast<int>(std::floor(position(dim) / rc));
        nidx[dim] = std::min(nidx[dim], nmax[dim]-1);
        nidx[dim] = std::max(nidx[dim], 0);
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
   * TODO include a distinction for the cutoff: with respect to the i-atom only
   * or with respect to the j,k,l etc. atom. I.e. the cutoff goes with the
   * Order?
   */
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
    n_i_atoms = this->manager.get_size();
    n_j_atoms = 0;
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

    // using Manager_t = StructureManager<ManagerImplementation>;
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
          auto j_add = pair.back();
          if (j_add > i_atoms.back()) {
            current_j_atoms.insert(j_add);
          }
        }
      }

      //! delete existing cluster atoms from list to build additional neighbours
      std::vector<size_t> atoms_to_add{};
      std::set_difference(current_j_atoms.begin(), current_j_atoms.end(),
                          current_i_atoms.begin(), current_i_atoms.end(),
                          std::inserter(atoms_to_add, atoms_to_add.begin()));

      manager.add_entry_number_of_neighbours();
      if (atoms_to_add.size() > 0) {
        for (auto j: atoms_to_add) {
          manager.add_neighbour_of_cluster(j);
        }
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Builds full neighbour list. Triclinicity is accounted for. The general idea
   * is anchor a mesh the origin of the supplied cell. Then the mesh is extended
   * into space until it is as big as the maximum cell coordinate plus one
   * cutoff. This mesh has boxes of size cutoff. Depending on the periodicity of
   * the mesh, ghost atoms are added by shifting all i-atoms by the cell
   * vectors. All i-atoms and the ghost atoms are then sorted into the
   * respective boxes of the mesh and a stencil anchored at the boxes of the
   * i-atoms is used to build the atom neighbourhoods from the 9 (2d) or 27 (3d)
   * boxes, which have to be checked for the neighbourhood . Correct periodicity
   * is ensured by the placement of the ghost atoms. The resulting neighbourlist
   * is full and not strict.
   */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::make_full_neighbour_list() {
    using Vector_t = Eigen::Matrix<double, traits::Dim, 1>;

    //! short hands for variable
    size_t natoms{this->manager.size()};
    const auto dim{traits::Dim};
    auto periodicity = this->manager.get_periodic_boundary_conditions();

    std::cout << periodicity[0] << " "
              << periodicity[1] << " "
              << periodicity[2] << " " << std::endl;

    auto cell{manager.get_cell()};
    double cutoff{this->cutoff};
    std::array<int, dim> nboxes_per_dim;

    //! vector for storing the atom indices of each box
    std::vector<std::vector<int>> atoms_in_box{};

    /**
     * minimum/maximum coordinate of mesh for neighbour list; depends on
     * cell triclinicity and cutoff
     */
    Vector_t mesh_min(dim);
    Vector_t mesh_max(dim);
    Vector_t cell_max(dim);
    mesh_min.setZero();
    mesh_max.setZero();

    //! maximum cell position
    for (auto i{0}; i < dim; ++i) { cell_max[i] = cell.col(i)[i];}

    //! calculate minimum mesh coordinates
    for (auto i{0}; i < dim; ++i) {
      mesh_min[i] -= cutoff;
    }

    //! calculate maximum mesh coordinates
    for (auto i{0}; i < dim; ++i) {
      auto dx =  cell.col(i)[i] + 2. * cutoff;
      std::cout << "dx " << dx << std::endl;
      int n(std::ceil(dx / cutoff));
      //! max is mesh origin + dx
      mesh_max[i] = mesh_min[i] + n * cutoff;
      nboxes_per_dim[i] = n;
    }

    /**
     * TODO possible future optimization for cells large triclinicity: use
     * triclinic coordinates and explicitly check for the 'skin' around the cell
     */
    //! generate ghost atom indices and position
    for (auto atom : this->get_manager()) {
      auto pos = atom.get_position();
      for (auto i{0}; i < dim; ++i) {
        if(periodicity[i]) {
          //! exclude m=0, because it is the i atom itself
          for(auto m{-1}; m < 2; m+=2) {
            //! shift i atom along cell vector i
            auto pos_ghost = pos + cell.col(i)*m;

            //! shift position to mesh origin
            auto pos_lower  = pos_ghost.array() - mesh_min.array();
            auto pos_greater  = pos_ghost.array() - mesh_max.array();
            //! shifted position inside mesh?
            auto f_lt = (pos_lower.array() > 0.).all();
            auto f_gt = (pos_greater.array() < 0.).all();
            if (f_lt and f_gt) {
              //! nex atom index is size, since start is at index = 0
              auto new_atom_index = this->get_size_with_ghosts();
              // new_atom_index++;
              this->add_ghost_atom(new_atom_index, pos_ghost);
            }
          }
        }
      }
    }
    std::cout << "natoms " << natoms
              << " nghosts " << n_j_atoms << std::endl;

    //! neighbour boxes
    internal::IndexContainer<dim> atom_id_cell{nboxes_per_dim};

    // i-atoms sorting into boxes
    for (size_t i{0}; i < this->n_i_atoms; ++i) {
      auto pos = this->get_position(i);
      auto dpos = pos - mesh_min;
      auto idx = internal::get_box_index(dpos, cutoff, nboxes_per_dim);
      atom_id_cell[idx].push_back(i);
    }

    //! ghost atoms sorting into boxes
    for (size_t i{0}; i < this->n_j_atoms; ++i) {
      auto ghost_pos = this->get_ghost_position(i);
      auto dpos = ghost_pos - mesh_min;
      auto idx  = internal::get_box_index(dpos, cutoff, nboxes_per_dim);
      auto ghost_atom_index = i + this->n_i_atoms;
      atom_id_cell[idx].push_back(ghost_atom_index);
    }

    //! go through atoms and build neighbour list
    int offset{0};
    for (size_t i{0}; i < this->n_i_atoms; ++i) {
      int nneigh{0};
      auto pos = this->get_position(i);
      auto dpos = pos - mesh_min;
      auto idx = internal::get_box_index(dpos, cutoff, nboxes_per_dim);
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
  //! Extend an existing neighbour list to the next order
  template <class ManagerImplementation>
  template <size_t MaxOrder, bool IsDummy>
  struct AdaptorMaxOrder<ManagerImplementation>::IncreaseMaxOrder {

    static void increase_maxorder(AdaptorMaxOrder<ManagerImplementation>
                                  * manager_max) {

      static_assert(MaxOrder > 2,
                    "no neighbourlist present; extension not possible.");

      for (auto atom : manager_max->manager) {
        //! Order 1, Order variable is at 0, atoms, index 0
        using AddOrderLoop = AddOrderLoop<atom.order(),
                                          atom.order() == traits::MaxOrder-1>;
        auto & atom_cluster_indices{std::get<0>
            (manager_max->cluster_indices_container)};
        atom_cluster_indices.push_back(atom.get_cluster_indices());
        AddOrderLoop::loop(atom, *manager_max);
      }

      //! correct the offsets for the new cluster order
      manager_max->set_offsets();
      //! add correct cluster_indices for the highest order
      auto & max_cluster_indices
      {std::get<traits::MaxOrder-1>(manager_max->cluster_indices_container)};
      max_cluster_indices.fill_sequence();
    }
  };

  template <class ManagerImplementation>
  template <size_t MaxOrder>
  struct AdaptorMaxOrder<ManagerImplementation>::
  IncreaseMaxOrder<MaxOrder, true> {
    static void increase_maxorder(AdaptorMaxOrder<ManagerImplementation>
                                  * manager_max) {
      manager_max->nb_neigh.resize(0);
      manager_max->offsets.resize(0);
      manager_max->neighbours.resize(0);
      manager_max->make_full_neighbour_list();
    }
  };


  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::update() {
    /**
     * Standard case, increase an existing neighbour list or triplet list to a
     * higher Order
     */
    using IncreaseMaxOrder = IncreaseMaxOrder<traits::MaxOrder,
                                              (traits::MaxOrder==2)>;
    IncreaseMaxOrder::increase_maxorder(this);
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
  inline size_t AdaptorMaxOrder<ManagerImplementation>::
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
    } else if (Order == traits::MaxOrder-1) {
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
       * If not accessible at this order, call lower Order offsets from lower
       * order manager(s).
       */
      return this->manager.get_offset_impl(counters);
    }
  }
}  // rascal

#endif /* ADAPTOR_MAXORDER_H */

/**
 * TODO: The construction of triplets is fine, but they occur multiple times. We
 * probably need to check for increasing atomic index to get rid of
 * duplicates. But this is in general a design decision, if we want full/half
 * neighbour list and full/half/whatever triplets and quadruplets
 */
