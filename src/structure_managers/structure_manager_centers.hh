/**
 * file   structure_manager_centers.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   06 August 2018
 *
 * @brief Manager with atoms and centers
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef STRUCTURE_MANAGER_CENTERS_H
#define STRUCTURE_MANAGER_CENTERS_H

#include "structure_managers/structure_manager.hh"
#include "lattice.hh"
#include "basic_types.hh"
//! Some data types and operations are based on the Eigen library
#include <Eigen/Dense>

//! And standard header inclusion
#include <stdexcept>
#include <vector>
#include <array>

/**
 * All functions and classes are in the namespace <code>rascal</code>, which
 * ensures that they don't clash with other libraries one might use in
 * conjunction.
 */
namespace rascal {
  //! forward declaration for traits
  class StructureManagerCenters;

  //! traits specialisation for Centers manager

  /**
   * The traits are used for vector allocation and further down the processing
   * chain to determine what functionality the given StructureManager
   * already contains to avoid recomputation.  See also the implementation of
   * adaptors.
   */
  template <>
  struct StructureManager_traits<StructureManagerCenters> {
    constexpr static int Dim{3};
    constexpr static size_t MaxOrder{1}; //
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDirectionVectors{false};
    constexpr static bool HasDistances{false};
    using LayerByDimension = std::integer_sequence<size_t, 0>;
  };

  /**
   * Definition of the new StructureManager class. To add your own, please
   * stick to the convention of using 'StructureManagerYours', where 'Yours'
   * will give a hint of what it is about.
   */
  class StructureManagerCenters:
    // It inherits publicly everything from the base class
    public StructureManager<StructureManagerCenters>
  {
    /**
     * Publicly accessible variables and function of the class are given
     * here. These provide the interface to access the neighbourhood.
     */
  public:
    //! For convenience, the names are shortened
    using traits = StructureManager_traits<StructureManagerCenters>;
    using Parent = StructureManager<StructureManagerCenters>;
    //! Here you see why -- definition of used function return types
    using Vector_ref = typename Parent::Vector_ref;
    using AtomRef_t = typename Parent::AtomRef;

    /**
     *  Eigen::Map is a convenient way to access data in the 'Eigen-way', if it
     * is already stored in a contiguous array.  The positions of the JSON file
     * and the cell vectors are put into contiguous arrays, which are member
     * variables of this class. Access is provided via the Eigen::Maps
     *
     * The following types are defined to access the data. Since they cost
     * almost nothing to build, they are created on the fly via e.g. the
     * .get_positions() member function, if needed. Access to the cell vectors,
     * defined in the JSON file.
     */
    using Cell_t = Eigen::Matrix<double, traits::Dim,traits::Dim,Eigen::ColMajor>;
    using Cell_ref = Eigen::Map<Cell_t>;
    using AtomTypes_t = Eigen::Matrix<int, 1, Eigen::Dynamic>;
    using AtomTypes_ref = Eigen::Map<AtomTypes_t>;
    using PBC_t = Eigen::Matrix<bool, 1, traits::Dim>;
    using PBC_ref = Eigen::Map<PBC_t>;
    using Positions_t = Eigen::Matrix<double, traits::Dim,Eigen::Dynamic,Eigen::ColMajor>;
    using Positions_ref = Eigen::Map<Positions_t>;

    /**
     * Here, the types for internal data structures are defined, based on
     * standard types.  In general, we try to use as many standard types, where
     * possible. It reduces the dependance on external libraries. If you want to
     * use e.g. the <code>Eigen</code> library, the access to the data can be
     * given via the above mentioned Eigen::Maps, which wrap arround contiguous
     * arrays of the internally saved data. It should be straight forward to
     * even include native <code>Eigen</code> types, since the compilation
     * checks for the library by default.
     */

    /**
     * A ClusterRef_t is a return type for iterators. It gives a light-weight
     * reference to an atom, a pair, a triplet,... to the AtomRefs of all
     * implicated atoms.  The template parameters Order and MaxOrder give the
     * pair/triplet/ and the maximum body order, e.g. up to pair level.
     * To increase the MaxOrder, use an <code>adaptor</code>.
     */
    template <size_t Order>
    using ClusterRef_t = typename Parent::template ClusterRef<Order>;

    //! Default constructor
    StructureManagerCenters() = default;

    //! Copy constructor
    StructureManagerCenters(const StructureManagerCenters &other) = delete;

    //! Move constructor
    StructureManagerCenters(StructureManagerCenters &&other) = default;

    //! Destructor
    virtual ~StructureManagerCenters() = default;

    //! Copy assignment operator
    StructureManagerCenters&
    operator=(const StructureManagerCenters &other) = delete;

    //! Move assignment operator
    StructureManagerCenters&
    operator=(StructureManagerCenters &&other) = default;

    /**
     * This member function invokes the reinitialisation of data. E.g. when the
     * atom positions are provided by a simulation method, which evolves in
     * time, this function updates the data. In this example, .update() takes no
     * arguments, because it relies on the data provided by a file. It is read
     * by invoking .read_structure_from_json(). Invoking update also builds a
     * full and half neighbour list.  The update function is required, every
     * time the list changes. It is implemented here with a dependency to the
     * JSON interface. An update of the positions in the JSON object
     * <code>atoms_object</code> needs to be followed by a call to this
     * function.
     * @param cutoff Property, which defines the cutoff in the
     * neighbourlist. Can in the future be combined with cutoff_skin for a
     * Verlet type list
     */
    void update(const Eigen::Ref<const Eigen::MatrixXd> positions,
                const Eigen::Ref<const VecXi>  particle_types,
                const Eigen::Ref<const VecXi> center_ids,
                const Eigen::Ref<const Eigen::MatrixXd> cell,
                const std::array<bool,3>& pbc);

    void build(const Eigen::Ref<const Eigen::MatrixXd> positions,
                const Eigen::Ref<const VecXi>  particle_types,
                const Eigen::Ref<const VecXi> center_ids,
                const Eigen::Ref<const Eigen::MatrixXd> cell,
                const std::array<bool,3>& pbc);

    //1 required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    // void reset_impl(const int & natoms);
    // // TODO

    /**
     * Returns a traits::Dim by traits::Dim matrix with the cell vectors of the
     * structure.
     */
    inline Cell_ref get_cell() {
      return Cell_ref(this->cell.data(), traits::Dim,
                      this->cell.size()/traits::Dim);
    }

    //! Returns the type of a given atom, given an AtomRef
    inline int get_atom_type(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto t = this->get_atom_types();
      return t(index);
    }

    //! Returns an a map with all atom types.
    inline AtomTypes_ref get_atom_types() {
      return AtomTypes_ref(this->particle_types.data(),
                           this->particle_types.size());
    }

    //! Returns a map of size traits::Dim with 0/1 for periodicity
    inline PBC_ref get_periodic_boundary_conditions() {
      return PBC_ref(this->pbc.data());
    }

    //! Returns the position of an atom, given an AtomRef
    inline Vector_ref get_position(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto p = this->get_positions();
      auto * xval{p.col(index).data()};
      return Vector_ref(xval);
    }

    //! Returns the position of an atom, given an atom index
    inline Vector_ref get_position(const size_t & atom_index) {
      auto p = this->get_positions();
      auto * xval{p.col(atom_index).data()};
      return Vector_ref(xval);
    }

    /**
     * Returns the position of a neighbour. In case of periodic boundary
     * conditions, the get_neighbour_position should return a different
     * position, if it is a ghost atom.
     */
    template<size_t Order, size_t Layer>
    inline void get_neighbour_position(const ClusterRefKey<Order, Layer>
                                             & ) {

      static_assert(true,
                    "this implementation only work with atoms.");

    }

    //! returns a map to all atomic positions.
    inline Positions_ref get_positions() {
      return Positions_ref(this->positions.data(), traits::Dim,
                           this->positions.size()/traits::Dim);
    }

    //! returns number of I atoms in the list
    inline size_t get_size() const {
      return this->natoms;
    }

    //! returns the number of neighbours of a given i atom
    template<size_t Order, size_t Layer>
    inline void get_cluster_size(const ClusterRefKey<Order, Layer>
                                   & ) const {
      static_assert(true,
                    "this implementation only handles atoms.");
    }

    //! Cluster size is the number of neighbours here
    inline size_t get_cluster_size(const int & ) const {
      return 1;
    }

    /**
     * Return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
                                  & counters) const;

    //! Function for returning the number of atoms, pairs, tuples, etc.
    size_t get_nb_clusters(size_t cluster_size) const;



  protected:

    /**
     * Since the data from the <code>atoms_object</code>, especially the
     * positions are not contiguous in memory (they are
     * <code>std:vector<std::vector<double>></code>) and those are the ones,
     * which are used heavily, they are put into contiguous memory
     * vector. Access to those are provided via Eigen:Map.
     */

    std::vector<AtomRef_t> particles;
    std::vector<AtomRef_t> centers; //!
    Positions_t positions; //!
    AtomTypes_t particle_types;
    Lattice lattice;
    Cell_t cell; // to simplify get_neighbour_position()
    std::array<bool,traits::Dim> pbc;
    /**
     * A vector which stores the absolute offsets for each atom to access the
     * correct variables in the neighbourlist.
     */
    std::vector<size_t> offsets{};

    /**
     * Convenience variables, which are set during <code>update</code> or while
     * the neighbourlist is built
     */
    //! Total number of atoms in structure
    size_t natoms{};


    /**
     * A switch to make the class verbose and give screen output about its
     * processes.
     */
    constexpr static bool verbose{false};

  private:
  };

  /* ---------------------------------------------------------------------- */
  // used for buildup
  template<size_t Order>
  inline size_t StructureManagerCenters::
  get_offset_impl(const std::array<size_t, Order> & counters) const {
    // TODO: Check this static_assert for validity
    // static_assert (Order == 1, "this manager can only give the offset "
    //                "(= starting index) for a pair iterator, given the i atom "
    //                "of the pair");
    return this->offsets[counters.front()];
  }

}  // rascal

#endif /* STRUCTURE_MANAGER_CENTERS_H */
