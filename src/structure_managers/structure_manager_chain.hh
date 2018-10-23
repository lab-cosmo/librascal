/**
 * file   structure_manager_chain.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   30 May 2018
 *
 * @brief Neighbourhood manager for polyalanine chain, reading
 *        structure from json file
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

/** NOTE to Developer
 * -----------------
 * This implementation of StructureManager is designed to be a
 * developer tutorial. It reads a molecular structure from a JSON file
 * in the format of ASE (Atomic Simulation environment,
 * https://wiki.fysik.dtu.dk/ase/index.html), transfers the data into
 * efficient containers where necessary and builds a full as well as a
 * half neighbour list based on the linked cell / linked list
 * algorithm.  Further processing of the structure is done via
 * adaptors, which are explained elsewhere.
 * The code is heavily commented to explain the thought process behind
 * the library design. Fell free to use this example as a basis for
 * your own implementation of a StructureManager.
 * Some very basic knowledge of C++ is needed, but the idea of this
 * tutorial is also to explain some recent features of the language
 * and why they are used here. Please also read the coding convention
 * to adhere to the coding style used throughout the library.
 */

//! Always use header guards
#ifndef STRUCTURE_MANAGER_CHAIN_H
#define STRUCTURE_MANAGER_CHAIN_H

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
  class StructureManagerChain;

  //! traits specialisation for Chain manager

  /**
   * The traits are used for vector allocation and further down the processing
   * chain to determine what functionality the given StructureManager
   * already contains to avoid recomputation.  See also the implementation of
   * adaptors.
   */
  template <>
  struct StructureManager_traits<StructureManagerChain> {
    constexpr static int Dim{3};
    constexpr static size_t MaxOrder{2}; //
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDirectionVectors{false};
    constexpr static bool HasDistances{false};
    using LayerByOrder = std::integer_sequence<size_t, 0, 0>;
  };

  /**
   * Definition of the new StructureManager class. To add your own, please
   * stick to the convention of using 'NeighbourhoofManagerYours', where 'Yours'
   * will give a hint of what it is about.
   */
  class StructureManagerChain:
    // It inherits publicly everything from the base class
    public StructureManager<StructureManagerChain>
  {
    /**
     * Publicly accessible variables and function of the class are given
     * here. These provide the interface to access the neighbourhood.
     */
  public:
    //! For convenience, the names are shortened
    using traits = StructureManager_traits<StructureManagerChain>;
    using Parent = StructureManager<StructureManagerChain>;
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


    // TODO change the ref to use the types defined in AtomicStructure
    using Cell_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
                                              Eigen::Dynamic>>;
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
    using NeighbourList_t = std::vector<std::vector<int>>;
    using HalfNeighbourList_t = std::vector<int>;
    using NumNeigh_t = std::vector<int>;
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
    StructureManagerChain() = default;

    //! Copy constructor
    StructureManagerChain(const StructureManagerChain &other) = delete;

    //! Move constructor
    StructureManagerChain(StructureManagerChain &&other) = default;

    //! Destructor
    virtual ~StructureManagerChain() = default;

    //! Copy assignment operator
    StructureManagerChain&
    operator=(const StructureManagerChain &other) = delete;

    //! Move assignment operator
    StructureManagerChain&
    operator=(StructureManagerChain &&other) = default;

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
    void update(double cutoff);

    //1 required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    // void reset_impl(const int & natoms);
    // // TODO

    /**
     * Returns a traits::Dim by traits::Dim matrix with the cell vectors of the
     * structure.
     */
    inline Cell_ref get_cell() {
      return Cell_ref(this->cell_data.data(), traits::Dim,
                      this->cell_data.size()/traits::Dim);
    }

    //! Returns the type of a given atom, given an AtomRef
    inline int get_atom_type(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto t = this->get_atom_types();
      return t(index);
    }

    //! Returns an a map with all atom types.
    inline AtomTypes_ref get_atom_types() {
      return AtomTypes_ref(this->atoms_object.type.data(),
                           this->atoms_object.type.size());
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

    /**
     * Returns the position of a neighbour. In case of periodic boundary
     * conditions, the get_neighbour_position should return a different
     * position, if it is a ghost atom.
     */
    template<size_t Order, size_t Layer>
    inline Vector_ref get_neighbour_position(const ClusterRefKey<Order, Layer>
                                             & cluster) {
      static_assert(Order > 1,
                    "Only possible for Order > 1.");
      static_assert(Order <= traits::MaxOrder,
                    "this implementation should only work up to MaxOrder.");
      return this->get_position(cluster.back());
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

    //! returns the number of neighbours of a given i atom
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
                                   & cluster) const {
      // TODO: Check for <= or < ?!
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles atoms and pairs.");
      return this->numneigh[cluster.back()];
    }

    //! return the index-th neighbour of cluster
    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
                                     & cluster,
                                     size_t index) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms and pairs.");
      auto && i_atom_id{cluster.back()};
      auto && off{this->offsets[i_atom_id]};
      return this->halfneigh[off + index];
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

    /**
     * A helper-function for the linked cell algorithm. It is used for the
     * division of the box into the cells.
     */
    inline double get_box_length(int dimension);

    /**
     * Function for constructing a full neighbourlist as well as a half
     * neighbourlist. Since the JSON file does not provide it and it is
     * necessary for calculating follow-up properties.
     */
    void make_neighbourlist(double cutoff);

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
    //! Number of pairs
    size_t nb_pairs{};
    //! A list of atom indeces (int), here it is just the number in the list
    Ilist_t ilist{};
    /**
     *  A full neighbour list. Each entry in allneigh is a vector of neighbours
     * of the given atom id
     */
    NeighbourList_t allneigh{};

    /**
     * Half neighbourlist. It is a contiguous vector where all neighbours are
     * listed once. Use in conjunction with <code>numneigh</code> to get the
     * number of neighbours per atom. This can be used to calculate a property
     * like the atomic distance.
     */
    HalfNeighbourList_t halfneigh{};
    /**
     * A vector which hold the number of neighbours (half neighbourlist) for
     * each atom.
     */
    NumNeigh_t numneigh{};

    //! Skin is an additional layer for the neighbour list (Verlet)
    double cutoff_skin{0.0};
    /**
     * A vector which stores the absolute offsets for each atom to access the
     * correct variables in the neighbourlist.
     */
    std::vector<size_t> offsets{};

    //! For linked cell algorithm; saved for possible later use.
    std::vector<int> ll{};
    std::vector<int> lc{};

    /**
     * Given a position and a discretization of the volume, this function gives
     * the traits:Dim box coordinates of the position in the cell. Usage for
     * example in construction of neighbour list
     */
    inline std::vector<int>
    get_box_index(Vector_ref& position,
                  std::vector<double>& rc,
                  Eigen::Matrix<double, 1, traits::Dim> offset,
                  std::vector<int> nmax);

    /**
     * Function which collects the neighbour atom id and writes it to the
     * member variable <code>allneigh</code>.
     */
    inline void collect_neighbour_info_of_atom(const int i,
                                               const std::vector<int> boxidx,
                                               const std::vector<int> nmax);

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
  inline size_t StructureManagerChain::
  get_offset_impl(const std::array<size_t, Order> & counters) const {
    /**
     * The static assert with <= is necessary, because the template parameter
     * ``Order`` is one Order higher than the MaxOrder at the current
     * level. The return type of this function is used to build the next Order
     * iteration.
     */
    static_assert (Order <= traits::MaxOrder, "this manager can only give the"
                   " offset (= starting index) for a pair iterator, given the"
                   " i atom of the pair");
    return this->offsets[counters.front()];
  }

}  // rascal

#endif /* STRUCTURE_MANAGER_CHAIN_H */
