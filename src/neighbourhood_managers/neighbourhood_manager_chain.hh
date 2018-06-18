/**
 * file   neighbourhood_manager_chain.hh
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

/* NOTE to Developer
 * -----------------
 * This implementation of NeighbourhoodManager is designed to be a
 * developer tutorial. It reads a molecular structure from a JSON file
 * in the format of ASE (Atomic Simulation environment,
 * https://wiki.fysik.dtu.dk/ase/index.html), transfers the data into
 * efficient containers where necessary and builds a full as well as a
 * half neighbour list based on the linked cell / linked list
 * algorithm.  Further processing of the structure is done via
 * adaptors, which are explained elsewhere.
 * The code is heavily commented to explain the thought process behind
 * the library design. Fell free to use this example as a basis for
 * your own implementation of a NeighbourhoodManager.
 * Some very basic knowledge of C++ is needed, but the idea of this
 * tutorial is also to explain some recent features of the language
 * and why they are used here. Please also read the coding convention
 * to adhere to the coding style used throughout the library.
 */

// Always use header guards
#ifndef NEIGHBOURHOOD_MANAGER_CHAIN_H
#define NEIGHBOURHOOD_MANAGER_CHAIN_H

// Each actual implementation of a NeighbourhoodManager is based
// on the given interface
#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/property.hh"

// An external header-library/header-class, which makes it easy to use
// the JSON as a first class data type. See
// https://github.com/nlohmann/json for documentation.
#include "json.hpp"

// Some data types and operations are based on the Eigen library
#include <Eigen/Dense>

// And standard header inclusion
#include <stdexcept>
#include <vector>

// For convenience
using json = nlohmann::json;

// All functions and classes are in the namespace <code>rascal</code>,
// which ensures that they don't clash with other libraries one might
// use in conjunction.
namespace rascal {

  //
  namespace JSONTransfer {

    // To read from a JSON file and deserialize the content, the used
    // class needs a <code>struct</code> with standard data types.
    struct AtomicStructure {
      /**
	 \param cell is a vector a vector of vectors which holds the cell unit
	 vectors.
	 \param type a vector of integers which holds the atomic type
	 (coordination number).
	 \param pbc is a 0/1 vector which says, where periodic boundary
	 conditions are applied.
	 \param position is a vector of vectors which holds the atomic
	 positions.
       */
      std::vector<std::vector<double>> cell{};
      std::vector<int> type{};
      std::vector<int> pbc{};
      std::vector<std::vector<double>> position{};
    };

    // This function is used to convert to the JSON format with the
    // given keywords. It is an overload of the function defined in the header
    // class json.hpp.
    // Inline needed, otherwise it is a multiple definition
    inline void to_json(json & j, AtomicStructure& s) {
      j = json{
	{"cell", s.cell},
	{"numbers", s.type},
	{"pbc", s.pbc},
	{"positions", s.position}
      };
    }

    // This function is used to read from the JSON file and convert
    // the data into standard types. It is an overload of the function defined
    // in json.hpp class header.
    inline void from_json(const json& j, AtomicStructure& s) {
      s.cell = j.at("cell").get<std::vector<std::vector<double>>>();
      s.type = j.at("numbers").get<std::vector<int>>();
      s.pbc = j.at("pbc").get<std::vector<int>>();
      s.position = j.at("positions").get<std::vector<std::vector<double>>>();
    }
  }

  //! forward declaration for traits
  class NeighbourhoodManagerChain;

  //! traits specialisation for Chain manager

  // The traits are used for vector allocation and further down the
  // processing chain to determine what the given Neighbourhoodmanager
  // already contains to avoid recomputation. See also the
  // implementation of adaptors.
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerChain> {
    constexpr static int Dim{3};
    constexpr static int MaxLevel{2}; //
  };
  // Definition of the new NeighbourhoodManager class. To add your
  // own, please stick to the convention of using
  // 'NeighbourhoofManagerYours', where 'Yours' will give a hin of
  // what it is about.
  class NeighbourhoodManagerChain:
    // It inherits publicly everything from the base class
    public NeighbourhoodManagerBase<NeighbourhoodManagerChain>
  {
    // Publicly accessible variables and function of the class are
    // given here. These provide the interface to access the neighbourhood.
  public:
    // For convenience, the names are shortened
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerChain>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerChain>;
    // Here you see why -- definition of used function return types
    using Vector_ref = typename Parent::Vector_ref;
    using AtomRef_t = typename Parent::AtomRef;

    // Eigen::Map is a convenient way to access data in the
    // 'Eigen-way', if it is already stored in a contiguous array.
    // The positions of the JSON file and the cell vectors are put
    // into contiguous arrays, which are member variables of this
    // class. Access is provided via the Eigen::Maps

    // The following types are defined to access the data. Since they
    // cost almost nothing to build, they are created on the fly via
    // e.g. the .get_positions() member function, if needed.
    // Access to the cell vectors, defined in the JSON file.
    using Cell_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
					      Eigen::Dynamic>>;
    // Access to the atom types, defined in the JSON file.
    using AtomTypes_ref = Eigen::Map<Eigen::Matrix<int, 1, Eigen::Dynamic>>;
    // Access to periodic boundary conditions, defined in JSON file, x,y,z
    using PBC_ref = Eigen::Map<Eigen::Matrix<int, 1, traits::Dim>>;
    // Access to an array of all given positions
    using Positions_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
						  Eigen::Dynamic>>;

    // Here, the types for internal data structures are defined, based on
    // standard types.  In general, we try to use as many standard types, where
    // possible. It reduces the dependance on external libraries. If you want to
    // use e.g. the <code>Eigen</code> library, the access to the data can be
    // given via the above mentioned Eigen::Maps, which wrap arround contiguous
    // arrays of the internally saved data. It should be straight forward to
    // even include native <code>Eigen</code> types, since the compilation
    // checks for the library by default.
    using NeighbourList_t = std::vector<std::vector<int>>;
    using HalfNeighbourList_t = std::vector<int>;
    using NumNeigh_t = std::vector<int>;
    using Ilist_t = std::vector<int>;
    // A ClusterRef_t is a return type for iterators. It gives a
    // light-weight reference to an atom, a pair, a triplet,... to the
    // AtomRefs of all implicated atoms.
    // The template parameters Level and MaxLevel give the
    // pair/triplet/ and the maximum depth, e.g. up to pair level.
    // To increase the MaxLevel, use an <code>adaptor</code>.
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template ClusterRef<Level, MaxLevel>;

    //! Default constructor
    NeighbourhoodManagerChain() = default;

    //! Copy constructor
    NeighbourhoodManagerChain(const NeighbourhoodManagerChain &other) = delete;

    //! Move constructor
    NeighbourhoodManagerChain(NeighbourhoodManagerChain &&other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerChain() = default;

    //! Copy assignment operator
    NeighbourhoodManagerChain&
    operator=(const NeighbourhoodManagerChain &other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerChain&
    operator=(NeighbourhoodManagerChain &&other) = default;

    // This member function invokes the reinitialisation of
    // data. E.g. when the atom positions are provided by a simulation
    // method, which evolves in time, this function updates the
    // data. In this example, .update() takes no arguments, because it
    // relies on the data provided by a file. It is read by invoking
    // .read_structure_from_json(). Invoking update also builds a full
    // and half neighbour list.
    /**
     * The update function is required, every time the list changes. It is
     * implemented here with a dependency to the JSON interface. An update of
     * the positions in the JSON object <code>atoms_object</code> needs to be
     * followed by a call to this function.
     * @param cutoff Property, which defines the cutoff in the
     *        neighbourlist. Can in the future be combined with cutoff_skin for a
     *        Verlet type list
     */
    void update(double cutoff);

    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    // void reset_impl(const int & natoms);
    // // TODO

    // Returns a traits::Dim by traits::Dim matrix with the cell
    // vectors of the structure.
    inline Cell_ref get_cell() {
      return Cell_ref(this->cell_data.data(), traits::Dim,
		      this->cell_data.size()/traits::Dim);
    }

    // Returns the type of a given atom.
    inline int get_atom_type(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto t = this->get_atom_types();
      return t(index);
    }

    // Returns an a map with all atom types.
    inline AtomTypes_ref get_atom_types() {
      return AtomTypes_ref(this->atoms_object.type.data(),
			   this->atoms_object.type.size());
    }

    // Returns a map of size traits::Dim with 0/1 for periodicity
    inline PBC_ref get_periodic_boundary_conditions() {
      return PBC_ref(this->atoms_object.pbc.data());
    }

    // Returns the position of an atom
    inline Vector_ref get_position(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto p = this->get_positions();
      auto * xval{p.col(index).data()};
      return Vector_ref(xval);
    }

    // Returns the position of a neighbour. In case of periodic
    // boundary conditions, the get_neighbour should return a
    // different position, if it is a ghost atom.
    template<int Level, int MaxLevel>
    inline Vector_ref get_neighbour_position(const ClusterRef_t<Level,
					     MaxLevel>& cluster) {
      static_assert(Level > 1,
		    "this implementation should only work with a neighbour");
      return this->get_position(cluster.get_atoms().back());
    }

    // Returns a map to all positions.
    inline Positions_ref get_positions() {
      return Positions_ref(this->pos_data.data(), traits::Dim,
			   this->pos_data.size()/traits::Dim);
    }

    // return number of I atoms in the list
    inline size_t get_size() const {
      return this->natoms;
    }

    // Returns the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level,
				   MaxLevel>& cluster) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs.");
      return this->numneigh[cluster.get_atoms().back().get_index()];
    }

    // Return the number of atoms forming the next higher cluster with
    // this one
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs.");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      return this->firstneigh[std::move(i_atom_id)][j_atom_id];
      return 0;
    }

    // Return unique id of a given atom. Here the atoms are numbered
    // from 0 to natoms (see .update() function).
    inline size_t get_atom_id(const Parent& /*cluster*/,
                              int i_atom_id) const {
      return this->ilist[i_atom_id];
    }

    /**
     * Return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<int Level, int MaxLevel>
    inline int
    get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    // Function for returning the number of atoms, pairs, tuples, etc.
    size_t get_nb_clusters(int cluster_size);

    // Function to read from a JSON file. Based on the above mentioned
    // header class
    void read_structure_from_json(const std::string filename);

    // A helper-function for the linked cell algorithm. It is used for
    // the division of the box into the cells.
    inline double get_box_length(int dimension);

    // Function for constructing a full neighbourlist as well as a
    // half neighbourlist. Since the JSON file does not provide it and
    // it is necessary for calculating follow-up properties.
    void make_neighbourlist(double cutoff);

    // Behind the interface.
  protected:

    // The variable for the JSON object is called
    // <code>atoms_object</code>, because ASE (see above) calls it's
    // read/write function an iterator for 'Atoms objects'.  It is
    // first class C++ data structure, which 'feels' like JSON.
    JSONTransfer::AtomicStructure atoms_object{};

    // Since the data from the <code>atoms_object</code>, especially
    // the positions are not contiguous in memory (they are
    // <code>std:vector<std::vector<double>></code>) and those are the
    // ones, which are used heavily, they are put into contiguous
    // memory vector. Access to those are provided via Eigen:Map.
    std::vector<double> cell_data{}; // volume cell dim*dim
    std::vector<double> pos_data{}; // array of positions dim*natoms

    // Convenience variables, which are set during <code>update</code>
    // or while the neighbourlist is built
    size_t natoms{}; // Total number of atoms in structure
    size_t nb_pairs{}; // Number of pairs
    Ilist_t ilist{}; // A list of atom indeces (int), here it is just
		   // the number in the list
    NeighbourList_t firstneigh{}; // A full neighbour list. Each entry
				  // in firstneigh is a vector of
				  // neighbours of the given atom id

    HalfNeighbourList_t halfneigh{}; // Half neighbourlist. It is a
				     // contiguous vector where all
				     // neighbours are listed
				     // once. Use in conjunction with
				     // <code>numneigh</code> to get
				     // the number of neighbours per
				     // atom. This can be used to
				     // calculate a property like the
				     // atomic distance.
    NumNeigh_t numneigh{}; // A vector which hold the number of
			   // neighbours (half neighbourlist) for each
			   // atom.

    double cutoff_skin{0.0}; // Added now, but intended for use later
			      // in a Verlet list.

    std::vector<int> offsets{}; // A vector which stores the absolute
				// offsets for each atom to access the
				// correct variables in the
				// neighbourlist.

    // For linked cell algorithm; saved for possible later use.
    std::vector<int> ll{}; // Linked list
    std::vector<int> lc{}; // Linear indexed linked cells

    // Given a position and a discretization of the volume, this
    // function gives the traits:Dim box coordinates of the position
    // in the cell. Usage for example in construction of neighbour
    // list
    inline std::vector<int>
    get_box_index(Vector_ref& position,
		  std::vector<double>& rc,
		  Eigen::Matrix<double, 1, traits::Dim> offset,
		  std::vector<int> nmax);

    // Function which collects the neighbour atom id and writes it to
    // the member variable <code>firstneigh</code>.
    inline void collect_neighbour_info_of_atom(const int i,
					       const std::vector<int> boxidx,
					       const std::vector<int> nmax);

    // A switch to make the class verbose and give screen output about
    // its processes.
    constexpr static bool verbose{false};

  private:
  };


  /* ---------------------------------------------------------------------- */
  // adjust for triplets
  template<int Level, int MaxLevel>
  inline int NeighbourhoodManagerChain::
  get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const {
    static_assert(Level == 2, "This class cas only handle single atoms and pairs");
    static_assert(MaxLevel == traits::MaxLevel, "Wrong maxlevel");

    auto atoms{cluster.get_atoms()};
    auto i{atoms.front().get_index()};
    auto j{cluster.get_index()};
    auto main_offset{this->offsets[i]};
    return main_offset + j;
  }

  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms
  template <>
  inline int NeighbourhoodManagerChain:: template
  get_offset_impl<1, 2>(const ClusterRef_t<1, 2>& cluster) const {
    return cluster.get_atoms().back().get_index();
  }


}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_CHAIN_H */
