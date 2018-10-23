/**
 * file   structure_manager_json.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   8 Aug 2018
 *
 * @brief structure manager for reading atomic structure from file
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

//! Always use header guards
#ifndef STRUCTURE_MANAGER_JSON_H
#define STRUCTURE_MANAGER_JSON_H

/**
 * Each actual implementation of a StructureManager is based on the given
 * interface ´structure_manager.hh´
 */
#include "structure_managers/structure_manager.hh"
#include "structure_managers/json_io.hh"

//! Some data types and operations are based on the Eigen library
#include <Eigen/Dense>

//! And standard header inclusion
#include <stdexcept>
#include <vector>

/**
 * All functions and classes are in the namespace <code>rascal</code>, which
 * ensures that they don't clash with other libraries one might use in
 * conjunction.
 */
namespace rascal {
  //! forward declaration for traits
  class StructureManagerJson;

  //! traits specialisation for Json manager

  /**
   * The traits are used for vector allocation and further down the processing
   * chain to determine what functionality the given StructureManager
   * already contains to avoid recomputation.  See also the implementation of
   * adaptors.
   */
  template <>
  struct StructureManager_traits<StructureManagerJson> {
    constexpr static int Dim{3};
    constexpr static size_t MaxOrder{1}; //
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDirectionVectors{false};
    constexpr static bool HasDistances{false};
    using LayerByOrder = std::integer_sequence<size_t, 0>;
  };

  /**
   * Definition of the new StructureManager class. To add your own, please
   * stick to the convention of using 'NeighbourhoofManagerYours', where 'Yours'
   * will give a hint of what it is about.
   */
  class StructureManagerJson:
    // It inherits publicly everything from the base class
    public StructureManager<StructureManagerJson>
  {
    /**
     * Publicly accessible variables and function of the class are given
     * here. These provide the interface to access the neighbourhood.
     */
  public:
    //! For convenience, the names are shortened
    using traits = StructureManager_traits<StructureManagerJson>;
    using Parent = StructureManager<StructureManagerJson>;
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
    using Cell_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
                                              traits::Dim>>; //Eigen::Dynamic>>;
    //! Access to the atom types, defined in the JSON file.
    using AtomTypes_ref = Eigen::Map<Eigen::Matrix<int, 1, Eigen::Dynamic>>;
    //! Access to periodic boundary conditions, defined in JSON file, x,y,z
    using PBC_ref = Eigen::Map<Eigen::Matrix<int, 1, traits::Dim>>;
    //! Access to an array of all given positions
    using Positions_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
                                                   Eigen::Dynamic>>;

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
    using Ilist_t = std::vector<int>;

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
    StructureManagerJson() = default;

    //! Copy constructor
    StructureManagerJson(const StructureManagerJson &other) = delete;

    //! Move constructor
    StructureManagerJson(StructureManagerJson &&other) = default;

    //! Destructor
    virtual ~StructureManagerJson() = default;

    //! Copy assignment operator
    StructureManagerJson&
    operator=(const StructureManagerJson &other) = delete;

    //! Move assignment operator
    StructureManagerJson&
    operator=(StructureManagerJson &&other) = default;

    /**
     * This member function invokes the reinitialisation of data. E.g. when the
     * atom positions are provided by a simulation method, which evolves in
     * time, this function updates the data. In this example, .update() takes no
     * arguments, because it relies on the data provided by a file. It is read
     * by invoking .read_structure_from_json().
     */
    void update();

    //! required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    // void reset_impl(const int & natoms);
    // // TODO

    /**
     * Returns a traits::Dim by traits::Dim matrix with the cell vectors of the
     * structure.
     */
    inline Cell_ref get_cell() {
      return Cell_ref(this->cell_data.data());
    }

    //! Returns the type of a given atom, given an AtomRef
    inline int get_atom_type(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto t = this->get_atom_types();
      return t(index);
    }

    //! Returns the type of a given atom, given an atom index
    inline int get_atom_type(const int & index) {
      auto t = this->get_atom_types();
      return t(index);
    }

    //! Returns an a map with all atom types.
    inline AtomTypes_ref get_atom_types() {
      return AtomTypes_ref(this->atoms_object.type.data(),
                           this->atoms_object.type.size());
    }

    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
                                   & /*cluster*/) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms.");
      return 1;
    }

    //! Returns a map of size traits::Dim with 0/1 for periodicity
    inline PBC_ref get_periodic_boundary_conditions() {
      return PBC_ref(this->atoms_object.pbc.data());
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

    //! returns a map to all atomic positions.
    inline Positions_ref get_positions() {
      return Positions_ref(this->pos_data.data(), traits::Dim,
                           this->pos_data.size()/traits::Dim);
    }

    //! returns number of I atoms in the list
    inline size_t get_size() const {
      return this->natoms;
    }

    //! return the index-th neighbour of cluster
    template<size_t Order, size_t Layer>
    // inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
    //                                  & cluster,
    //                                  size_t index) const {
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
                                     & ,
                                     size_t ) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms and pairs.");
      return 0;
    }

    //! return the atom_index of the index-th atom in manager
    inline int get_cluster_neighbour(const Parent& /*cluster*/,
                                     size_t index) const {
      return this->ilist[index];
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

    /**
     * Function to read from a JSON file. Based on the above mentioned header
     * class
     */
    void read_structure_from_json(const std::string filename);

    //! Behind the interface.
  protected:

    /**
     * The variable for the JSON object is called <code>atoms_object</code>,
     * because ASE (see above) calls it's read/write function an iterator for
     * 'Atoms objects'.  It is first class C++ data structure, which 'feels'
     * like JSON.
     */
    json_io::AtomicStructure atoms_object{};

    /**
     * Since the data from the <code>atoms_object</code>, especially the
     * positions are not contiguous in memory (they are
     * <code>std:vector<std::vector<double>></code>) and those are the ones,
     * which are used heavily, they are put into contiguous memory
     * vector. Access to those are provided via Eigen:Map.
     */
    std::vector<double> cell_data{};
    std::vector<double> pos_data{};

    /**
     * Convenience variables, which are set during <code>update</code> or while
     * the neighbourlist is built
     */
    //! Total number of atoms in structure
    size_t natoms{};
    //! A list of atom indeces (int), here it is just the number in the list
    Ilist_t ilist{};

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
  inline size_t StructureManagerJson::
  get_offset_impl(const std::array<size_t, Order> & counters) const {
    // TODO: Check this static_assert for validity
    static_assert (Order <= traits::MaxOrder, "this manager can not provide any"
                   " offsets.");
    return counters.front();
  }

}  // rascal

#endif /* STRUCTURE_MANAGER_JSON_H */
