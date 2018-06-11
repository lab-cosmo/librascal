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

#ifndef NEIGHBOURHOOD_MANAGER_CHAIN_H
#define NEIGHBOURHOOD_MANAGER_CHAIN_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/property.hh"

#include "json.hpp"

#include <Eigen/Dense>

#include <stdexcept>
#include <vector>
//#include <string>


using json = nlohmann::json;

namespace rascal {

  namespace JSONTransfer {

    //! JSON specific
    struct Molecule {
      std::vector<std::vector<double>> cell{};
      std::vector<int> type{};
      std::vector<bool> pbc{};
      std::vector<std::vector<double>> position{};
    };

    // resorted to inline, but not sure if that should be so?
    inline void to_json(json & j, Molecule& s) {
      j = json{
	{"cell", s.cell},
	{"numbers", s.type},
	{"pbc", s.pbc},
	{"positions", s.position}
      };
    }

    inline void from_json(const json& j, Molecule& s) {
      s.cell = j.at("cell").get<std::vector<std::vector<double>>>();
      s.type = j.at("numbers").get<std::vector<int>>();
      s.pbc = j.at("pbc").get<std::vector<bool>>();
      s.position = j.at("positions").get<std::vector<std::vector<double>>>();
    }
  }

  //! forward declaration for traits
  class NeighbourhoodManagerChain;

  //! traits specialisation for Chain manager
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerChain> {
    constexpr static int Dim{3};
    constexpr static int MaxLevel{3}; // triplets needed for angle
  };
  class NeighbourhoodManagerChain:
    public NeighbourhoodManagerBase<NeighbourhoodManagerChain>
  {
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerChain>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerChain>;
    using Vector_ref = typename Parent::Vector_ref;
    using AtomRef_t = typename Parent::AtomRef;

    // Maps for access to complete contiguous data sets as Eigen types
    /*
     * These are needed extra along with the accessors, because these
     * maps can not be member variables. But they cost 'nothing' to
     * to build and access to the complete dataset might be useful.
     */
    using Cell_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
					      Eigen::Dynamic>>;
    using AtomTypes_ref = Eigen::Map<Eigen::Matrix<int, 1, Eigen::Dynamic>>;
    using PBC_ref = Eigen::Map<Eigen::Matrix<int, 1, traits::Dim>>;
    using Positions_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
						  Eigen::Dynamic>>;

    using NeighbourList_t = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using NumNeigh_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;
    using Ilist_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;
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

    void update();
    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    //! something like this!
    // void reset_impl(const int & natoms);
    // // TODO
    // constexpr static int dim() {return traits::Dim;}


    // // return neighbour positions if Level > 1, temporary
    // inline Vector_ref get_neighbour_position(const AtomRef_t& atom,
    // 					     const AtomRef_t&,
    // 					     const int&) {
    //   // TODO
    // }

    // // return force vector

    inline Cell_ref get_cell() {
      return Cell_ref(this->cell_data.data(), traits::Dim,
		      this->cell_data.size()/traits::Dim);
    }

    inline int get_atom_type(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto t = this->get_atom_types();
      return t(index);
    }

    inline AtomTypes_ref get_atom_types() {
      return AtomTypes_ref(this->type_data.data(), this->type_data.size());
    }

    // TODO: invalid use of void expression?
    inline PBC_ref get_periodic_boundary_conditions() {
      return PBC_ref(this->pbc_data.data());
    }

    inline Vector_ref get_position(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto p = this->get_positions();
      auto * xval{p.col(index).data()};
      return Vector_ref(xval);
    }

    inline Positions_ref get_positions() {
      return Positions_ref(this->pos_data.data(), traits::Dim,
			   this->pos_data.size()/traits::Dim);
    }

    // return number of I atoms in the list
    inline size_t get_size() const {
      return this->natoms;
    }

    // return the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level,
				   MaxLevel>& cluster) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms, "
		    "pairs and triplets");
      return this->numneigh[cluster.get_atoms().back().get_index()];
    }

    // return the number of atoms forming the next higher cluster with this one
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const {
      // static_assert(Level == traits::MaxLevel-1,
      //               "this implementation only handles atoms, pairs "
      // 		    "and triplets");
      // auto && i_atom_id{cluster.get_atoms().back().get_index()};
      // return this->firstneigh[std::move(i_atom_id)][j_atom_id];
      return 0;
    }

    // return the number of neighbours of a given atom
    inline size_t get_atom_id(const Parent& /*cluster*/,
                              int i_atom_id) const {
      return this->ilist[i_atom_id];
    }

    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<int Level, int MaxLevel>
    inline int
    get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

    void read_structure_from_json(const std::string filename);

  protected:

    JSONTransfer::Molecule molecule_in;

    // References to contiguous vectors for nested types
    std::vector<double> cell_data{};
    std::vector<int> type_data{};
    // this is actually bool - but Eigen does not like it
    std::vector<int> pbc_data{};
    std::vector<double> pos_data{};

    // To be initialized by contruction of manager for actual neighbour use
    size_t natoms{}; // total number of atoms in structure
    // size_t nb_pairs{};
    // size_t nb_triplets{};
    // size_t nb_quadruplets{};
    Ilist_t ilist; // adhering to lammps-naming
    // adhering to lammps-naming: TODO will be initialized bei MaxLevel + 1 adaptor
    // NeighbourList_t firstneigh;
    NumNeigh_t numneigh; // adhering to lammps-naming

    std::vector<int> offsets;

  private:

  };


  /* ---------------------------------------------------------------------- */
  // adjust for triplets
  template<int Level, int MaxLevel>
  inline int NeighbourhoodManagerChain::
  get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const {
    static_assert(Level == 1,
		  "This class cas only handle single atoms; "
		  "use adaptors to increase MaxLevel.");
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
