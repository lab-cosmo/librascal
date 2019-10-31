/**
 * @file  src/structure_managers/structure_manager_centers.hh
 *
 * @ingroup group_structure_manager
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   06 August 2018
 *
 * @brief basic manager implementation with atoms and centers with ability to
 *        read from json file and be constructed from existing data
 *
 * Copyright  2018 Felix Musil, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_CENTERS_HH_
#define SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_CENTERS_HH_

// inclusion of librascal data structure, each manager is based on the interface
// given in `structure_manager.hh`
#include "atomic_structure.hh"
#include "basic_types.hh"
#include "json_io.hh"
#include "lattice.hh"
#include "structure_managers/structure_manager.hh"

// data types and operations are based on the Eigen library
#include <Eigen/Dense>

// standard header inclusion
#include <array>
#include <stdexcept>
#include <vector>

/**
 * All functions and classes are in the namespace <code>rascal</code>, which
 * ensures that they don't clash with other libraries one might use in
 * conjunction.
 */
namespace rascal {
  //! forward declaration for traits
  class StructureManagerCenters;

  /**
   * traits specialisation for ManagerCenters: traits are used for vector
   * allocation and further down the processing chain to determine what
   * functionality the given StructureManager already contains to avoid
   * recomputation. See also the implementation of adaptors.
   */
  template <>
  struct StructureManager_traits<StructureManagerCenters> {
    // this manager only works with 3D data (positions, cell)
    constexpr static int Dim{3};
    // MaxOrder is set to 1 because it has just atoms
    constexpr static size_t MaxOrder{1};
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDirectionVectors{false};
    constexpr static bool HasDistances{false};
    constexpr static bool HasCenterPair{false};
    constexpr static int StackLevel{0};
    using LayerByOrder = std::index_sequence<0>;
  };

  /**
   * StructureManagerCenters is an entry point to the neighbourlist. It takes
   * an atomic structure (positions, atomic number, cell and periodic boundary
   * conditions) and allows to select the atoms for which a neighbourlist will
   * be built (by default all atoms are considered). Note that all atoms are
   * considered as potential neighbors.
   *
   * The setting of the structure is done through the update function and with
   * the help of the AtomicStructure class.
   *
   * This manager allows for one level of iteration over the centers or i-atoms
   * that have been selected.
   */
  class StructureManagerCenters :
      // public inheritance of the base class
      public StructureManager<StructureManagerCenters>,
      public std::enable_shared_from_this<StructureManagerCenters> {
    // Publicly accessible variables and functions of the class are given
    // here. These provide the interface to access the structure and
    // subsequently calculated neighbourhood.
   public:
    // for convenience, the names are shortened
    using traits = StructureManager_traits<StructureManagerCenters>;
    using Parent = StructureManager<StructureManagerCenters>;
    // here you see why -- definition of used function return types
    using Vector_ref = typename Parent::Vector_ref;
    using AtomRef_t = typename Parent::AtomRef;
    using Children_t = typename Parent::Children_t;
    using ManagerImplementation_t = StructureManagerCenters;
    using ImplementationPtr_t = std::shared_ptr<StructureManagerCenters>;

    /**
     * Eigen::Map is a convenient way to access data in the 'Eigen-way', if it
     * is already stored in a contiguous array.  The positions of the JSON file
     * and the cell vectors are put into contiguous arrays, which are member
     * variables of this class. Access is provided via the Eigen::Maps
     *
     * The following types are defined globally in atomic_structure.hh as common
     * return types for accessing atom and cell related data.
     */
    using Cell_t = AtomicStructure<traits::Dim>::Cell_t;
    using Cell_ref = AtomicStructure<traits::Dim>::Cell_ref;

    using AtomTypes_t = AtomicStructure<traits::Dim>::AtomTypes_t;
    using AtomTypes_ref = AtomicStructure<traits::Dim>::AtomTypes_ref;
    using ConstAtomTypes_ref = AtomicStructure<traits::Dim>::ConstAtomTypes_ref;

    using PBC_t = AtomicStructure<traits::Dim>::PBC_t;
    using PBC_ref = AtomicStructure<traits::Dim>::PBC_ref;

    using Positions_t = AtomicStructure<traits::Dim>::Positions_t;
    using Positions_ref = AtomicStructure<traits::Dim>::Positions_ref;

    using ArrayB_t = AtomicStructure<traits::Dim>::ArrayB_t;
    using ArrayB_ref = AtomicStructure<traits::Dim>::ArrayB_ref;

    /**
     * Here, the types for internal data structures are defined, based on
     * standard types.  In general, we try to use as many standard types, where
     * possible. It reduces the dependance on external libraries. If you want to
     * use e.g. the <code>Eigen</code> library, the access to the data can be
     * given via the above mentioned Eigen::Maps, which wrap arround contiguous
     * arrays of the internally saved data. It should be straight forward to
     * even include native <code>Eigen</code> types, since the compilation
     * checks for the library by default.
     *
     * A ClusterRef_t is a return type for iterators. It gives a light-weight
     * reference to an atom, a pair, a triplet,... to the AtomRefs of all
     * implicated atoms.  The template parameters Order and MaxOrder give the
     * pair/triplet/ and the maximum body order, e.g. up to pair level.  To
     * increase the MaxOrder, use an <code>adaptor</code>.
     */
    template <size_t Order>
    using ClusterRef_t = typename Parent::template ClusterRef<Order>;

    //! Default constructor, values are set during .update() function
    StructureManagerCenters()
        : atoms_object{}, lattice{}, atoms_index{}, offsets{},
          n_center_atoms{}, natoms{} {}

    //! Copy constructor
    StructureManagerCenters(const StructureManagerCenters & other) = delete;

    //! Move constructor
    StructureManagerCenters(StructureManagerCenters && other) = delete;

    //! Destructor
    virtual ~StructureManagerCenters() = default;

    //! Copy assignment operator
    StructureManagerCenters &
    operator=(const StructureManagerCenters & other) = delete;

    //! Move assignment operator
    StructureManagerCenters &
    operator=(StructureManagerCenters && other) = delete;

    //! Updates the manager using the impl
    template <class... Args>
    void update(Args &&... arguments) {
      if (sizeof...(arguments) > 0) {
        // the structure has changed to tell it to the whole tree
        this->send_changed_structure_signal();
      }

      // update the underlying structure
      this->update_self(std::forward<Args>(arguments)...);
      this->set_update_status(true);

      // send the update signal to the tree
      this->update_children();
    }
    //! required for the construction of vectors, etc
    constexpr static int dim() { return traits::Dim; }

    /**
     * Returns a traits::Dim by traits::Dim matrix with the cell vectors of the
     * structure.
     */
    Cell_ref get_cell() { return Cell_ref(this->atoms_object.cell); }

    //! Returns the type of a given atom, given an AtomRef
    int get_atom_type(int atom_tag) const {
      auto && atom_index{this->get_atom_index(atom_tag)};
      return this->atoms_object.atom_types(atom_index);
    }

    //! Returns an a map with all atom types.
    ConstAtomTypes_ref get_atom_types() const {
      ConstAtomTypes_ref val(this->atoms_object.atom_types);
      return val;
    }

    //! Returns a map of size traits::Dim with 0/1 for periodicity
    PBC_ref get_periodic_boundary_conditions() {
      return PBC_ref(this->atoms_object.pbc);
    }

    ArrayB_ref get_center_atoms_mask() {
      return ArrayB_ref(this->atoms_object.center_atoms_mask);
    }

    //! Returns the position of an atom, given an AtomRef
    Vector_ref get_position(const AtomRef_t & atom) {
      auto atom_tag{atom.get_index()};
      return this->get_position(atom_tag);
    }

    //! Returns the position of an atom, given an atom tag
    Vector_ref get_position(int atom_tag) {
      auto && atom_index{this->get_atom_index(atom_tag)};
      auto p = this->get_positions();
      auto * xval{p.col(atom_index).data()};
      return Vector_ref(xval);
    }

    //! returns a map to all atomic positions.
    Positions_ref get_positions() {
      return Positions_ref(this->atoms_object.positions);
    }

    //! returns number of I atoms in the list
    size_t get_size() const { return this->n_center_atoms; }

    //! returns number of I atoms in the list, since at this level, center atoms
    //! and ghost atoms are not distinguishable.
    size_t get_size_with_ghosts() const { return this->n_center_atoms; }

    //! returns the number of neighbours of a given i atom
    template <size_t Order, size_t Layer>
    size_t get_cluster_size_impl(
        const ClusterRefKey<Order, Layer> & /*cluster*/) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms.");
      return 1;
    }

    template <size_t Order, size_t Layer>
    size_t get_atom_index(const ClusterRefKey<Order, Layer> & cluster) const {
      return this->get_atom_index(cluster.get_atom_tag());
    }

    size_t get_atom_index(int atom_tag) const {
      return std::get<0>(this->atoms_index)[atom_tag];
    }

    //! dummy function, since no neighbours are present her
    int get_neighbour_atom_tag(const Parent & /*parent*/, size_t index) const {
      // dummy argument is the atom itself, because if does not make sense at
      // this order
      return index;
    }

    //! Dummy function, since neighbours are not present at this Order
    template <size_t Order, size_t Layer>
    int get_neighbour_atom_tag(const ClusterRefKey<Order, Layer> & /*cluster*/,
                               size_t /*index*/) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms.");
      return 0;
    }

    int get_atom_tag(int raw_index) const { return raw_index; }

    size_t get_atom_tag(size_t raw_index) const { return raw_index; }

    /**
     * Return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template <size_t Order>
    size_t get_offset_impl(const std::array<size_t, Order> & counters) const;

    //! Function for returning the number of atoms
    size_t get_nb_clusters(size_t order) const;

    //! to follow the base class
    void update_self() {}

    /**
     * Use AtomObject to read the incoming structure
     */
    template <class... Args>
    void update_self(Args &&... arguments) {
      this->atoms_object.set_structure(std::forward<Args>(arguments)...);
      this->build();
    }

    const AtomicStructure<traits::Dim> & get_atomic_structure() const {
      return this->atoms_object;
    }

    size_t get_n_atoms() const { return this->natoms; }

   protected:
    //! makes atom tag lists and offsets
    void build();
    /**
     * Object which can interface to the json header to read and write atom
     * related data in the ASE format: positions, cell, periodicity, atom types
     * (corresponding to element numbers)
     */
    AtomicStructure<traits::Dim> atoms_object{};

    //! Lattice type for storing the cell and querying cell-related data
    Lattice<traits::Dim> lattice;

    /**
     * store atoms index per order,i.e.
     *   - atoms_index[0] lists all i-atoms that will be centered on and at the
     * back the atoms that won't be centered on are indexed too.
     *
     * The size of the manager is n_center_atoms so the loop over the center
     * atoms does not include the atoms that are not centered on.
     *
     * This array is expected to be accessed with atom_tags.
     *
     * It is important to index all atoms even if they are not centers since
     * they will potentially be neighbours and we want to be able to get the
     * atom_index of a neighbour.
     */
    std::array<std::vector<size_t>, traits::MaxOrder> atoms_index;

    /**
     * A vector which stores the absolute offsets for each atom to access the
     * correct variables in the neighbourlist. Here they are just the 1, 2, 3,
     * etc. corresponding to the sequence in which the atoms are provided
     * during construction.
     */
    std::vector<size_t> offsets{};

    //! Total number of center atoms in the structure
    size_t n_center_atoms{};

    //! Total number of atoms in the structure
    size_t natoms{};

    //! number of time the structure has been updated
    size_t n_update{0};
  };

  /* ---------------------------------------------------------------------- */
  //! used for construction (not for iteration)
  template <size_t Order>
  size_t StructureManagerCenters::get_offset_impl(
      const std::array<size_t, Order> & /*counters*/) const {
    static_assert(Order == 1, "This manager only handles atoms.");
    return 0;
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_CENTERS_HH_
