/**
 * file   neighbourhood_manager_Minimal.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager linked Minimal
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


#ifndef NEIGHBOURHOOD_MANAGER_MINIMAL_H
#define NEIGHBOURHOOD_MANAGER_MINIMAL_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/property.hh"
#include "lattice.hh"
#include "basic_types.hh"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

namespace rascal {

  //! forward declaration for traits
  class NeighbourhoodManagerMinimal;

  //! traits specialisation for Lammps manager
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerMinimal> {
    constexpr static int Dim{3};
    constexpr static int MaxLevel{1};
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
  };
  class NeighbourhoodManagerMinimal: public NeighbourhoodManagerBase<NeighbourhoodManagerMinimal>
  {
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerMinimal>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerMinimal>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using AtomRef_t = typename Parent::AtomRef;
    template <int Level>
    using ClusterRef_t = typename Parent::template ClusterRef<Level>;

    using AtomVectorField_t = Property<NeighbourhoodManagerMinimal, double, 1, 3>;

    //! Default constructor
    NeighbourhoodManagerMinimal()
    :particles{}, centers{} ,positions{},shifted_position{} ,lattice{}, cell{},pbc{} ,part2bin{} ,boxes{} ,number_of_neighbours{0} ,neighbour_bin_id{} , number_of_neighbours_stride{}, neighbour_atom_index{},particule_types{}
    {}

    //! Copy constructor
    NeighbourhoodManagerMinimal(const NeighbourhoodManagerMinimal &other) = delete;

    //! Move constructor
    NeighbourhoodManagerMinimal(NeighbourhoodManagerMinimal &&other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerMinimal() = default;

    //! Copy assignment operator
    NeighbourhoodManagerMinimal& operator=(const NeighbourhoodManagerMinimal &other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerMinimal& operator=(NeighbourhoodManagerMinimal &&other) = default;

    class Box;

    // return position vector
    inline Vector_ref get_position(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto * xval{this->positions.col(index).data()};
      return Vector_ref(xval);
    }

    inline Vector_ref get_shift(const int& i_bin_id, const int& shift_index);

    // return position vector
    // atom is the neighbour atom. center_atom is the current center. j_linear_id is the index of the current neighbour iterator.
    template<int Level>
    inline Vector_ref get_neighbour_position(const ClusterRef_t<Level>&
					     cluster) {
      static_assert(Level > 1,
		    "this implementation should only work with a neighbour");
      auto && j_linear_id = cluster.get_index();
      auto && i_atom_id{cluster.get_atoms().front().get_index()}; // center_atom index
      auto && i_bin_id{this->part2bin[i_atom_id]};
      auto && shift_index{this->neighbour_bin_id[i_bin_id][j_linear_id].get_index()};
      auto && j_atom_id{cluster.get_atoms().back().get_index()}; // neighbour atom index
      // TODO: find another way. This is a work around so that shifted_position lives longer than the function call but it is prone to side effects
      this->shifted_position = this->positions.col(j_atom_id) + this->cell * this->get_shift(i_bin_id,shift_index);
      auto * xval{this->shifted_position.col(0).data()};
      return Vector_ref(xval);
    }

    // return number of center in the list
    inline size_t get_size() const {
      return this->centers.size();
    }
    // return the id a given center atom
    inline size_t get_atom_id(const Parent& , int i_atom_id) const {
      return this->centers[i_atom_id].get_index();
    }

    //! return atom type
    inline int get_atom_type(const AtomRef_t& atom) {
      auto && index{atom.get_index()};
      return this->particule_types[index];
    }

    // return the index of the center corresponding to its neighbour image
    template<int Level>
    inline size_t get_atom_id(const ClusterRef_t<Level>& cluster, int j_atom_id) const {
      static_assert(Level <= traits::MaxLevel,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      auto && i_bin_id{this->part2bin[i_atom_id]};
      auto && ij_atom_id{this->neighbour_atom_index[i_bin_id][j_atom_id].get_index()};
      return ij_atom_id;
    }

    // return the number of neighbours of a given atom
    template<int Level>
    inline size_t get_cluster_size(const ClusterRef_t<Level>& cluster) const {
      static_assert(Level <= traits::MaxLevel,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      auto && box_id{this->part2bin[i_atom_id]};
      auto && size{this->neighbour_atom_index[box_id].size()};
      return size;
    }

    template<int Level>
    inline int get_offset_impl(const ClusterRef_t<Level>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

    void update(const Eigen::Ref<const Eigen::MatrixXd> positions,const Eigen::Ref<const VecXi>  particule_types,
               const Eigen::Ref<const VecXi> center_ids,
                const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);

    //Box get_box(const int& bin_id);

    size_t get_box_nb();

  protected:

    void build(const Eigen::Ref<const Eigen::MatrixXd> positions, const Eigen::Ref<const VecXi>  particule_types,
               const Eigen::Ref<const VecXi> center_ids,
               const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);

    void set_positions(const Eigen::Ref<const Eigen::MatrixXd> pos){
      this->positions = pos;
    }

    std::vector<AtomRef_t> particles;
    std::vector<AtomRef_t> centers; //!
    Matrix3XdC positions; //!
    Vector_t shifted_position;
    Lattice lattice;
    Cell_t cell; // to simplify get_neighbour_position()
    std::array<bool,3> pbc;
    std::vector<int> part2bin; //!
    std::vector<Box> boxes;
    size_t number_of_neighbours;
    std::vector<std::vector<AtomRef_t>> neighbour_bin_id;
    std::vector<size_t> number_of_neighbours_stride;
    std::vector<std::vector<AtomRef_t>> neighbour_atom_index;
    std::vector<int> particule_types;

  private:
  };

  /* ---------------------------------------------------------------------- */

  class NeighbourhoodManagerMinimal::Box {
  public:
    using Manager_t = NeighbourhoodManagerBase<NeighbourhoodManagerMinimal>;
    using AtomRef_t = typename NeighbourhoodManagerMinimal::AtomRef_t;
    using Vector_t = typename NeighbourhoodManagerMinimal::Vector_t;
    using Vector_ref = typename NeighbourhoodManagerMinimal::Vector_ref;
    //! Default constructor
    Box() = default;

    //! constructor
    Box(Manager_t& manager, const Vec3i_t& coord, const std::array<bool, 3>& pbc, const Vec3i_t& neigh_search, const Vec3i_t& nbins_c);
          //const std::array<std::array<Dim_t, 3>,2>& neigh_bounds,


    //! copy constructor
    Box(const Box & other) = default;
    //! assignment operator
    Box & operator=(const Box & other) = default;

    virtual ~Box() = default;

    constexpr static int dim() {return NeighbourhoodManagerMinimal::dim();}

    inline void push_particle_back(const int& part_index);

    inline size_t get_number_of_particles();

    inline size_t get_number_of_neighbours();

    inline size_t get_number_of_neighbour_box();

    inline void set_number_of_neighbours(const size_t& neigh_nb);

    inline  int get_neighbour_bin_index(const int& j_index);

    inline size_t get_particle_index(const int& index);

    inline Vector_ref get_neighbour_bin_shift(const int& neigh_bin_index){
      auto * xval{this->neighbour_bin_shift[neigh_bin_index].col(0).data()};
      return Vector_ref(xval);
    }

  protected:
    Manager_t & manager;
    std::vector<AtomRef_t> particles;
    std::vector<Vector_t,Eigen::aligned_allocator<Vector_t>> neighbour_bin_shift; //stores double for the dot product with the Cell vectors
    std::vector<int> neighbour_bin_index;
    size_t number_of_neighbours;
    Vec3i_t coordinates;
  };
  /* ---------------------------------------------------------------------- */

  template<int Level>
  inline int NeighbourhoodManagerMinimal::
  get_offset_impl(const ClusterRef_t<Level>& cluster) const {
    static_assert(Level == 2, "This class cas only handle single atoms and pairs");

    //auto atoms{};
    auto icenter{cluster.get_atoms().front().get_index()};
    auto stride{this->number_of_neighbours_stride[icenter]};
    auto j{cluster.get_index()};
    auto main_offset{stride+j};
    return main_offset;
  }

  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms

  template <>
  inline int NeighbourhoodManagerMinimal:: template
  get_offset_impl<1, 2>(const ClusterRef_t<1, 2>& cluster) const {
    return cluster.get_atoms().back().get_index();
  }
  /* ---------------------------------------------------------------------- */

  inline Vector_ref NeighbourhoodManagerMinimal::get_shift(const int& i_bin_id, const int& neigh_bin_index){
      return this->boxes[i_bin_id].get_neighbour_bin_shift(neigh_bin_index);
  }

  //----------------------------------------------------------------------------//
}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_MINIMAL_H */
