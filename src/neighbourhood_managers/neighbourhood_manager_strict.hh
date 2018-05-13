/**
 * file   neighbourhood_manager_cell.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager linked cell
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


#ifndef NEIGHBOURHOOD_MANAGER_STRICT_H
#define NEIGHBOURHOOD_MANAGER_STRICT_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/neighbourhood_manager_cell.hh"
#include "neighbourhood_managers/field.hh"
#include "basic_types.hh"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

namespace rascal {

  //! forward declaration for traits
  class NeighbourhoodManagerStrict;

  //! traits specialisation for Lammps manager
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerStrict> {
    constexpr static int Dim{3};
    constexpr static int MaxLevel{3};
  };
  class NeighbourhoodManagerStrict: public NeighbourhoodManagerBase<NeighbourhoodManagerStrict>
  {
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerStrict>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerStrict>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using AtomRef_t = typename Parent::AtomRef;
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template ClusterRef<Level, MaxLevel>;
    
    using AtomVectorProp_t = Field<NeighbourhoodManagerStrict, double, 3, 3>;
    using AtomScalarProp_t = Field<NeighbourhoodManagerStrict, double, 3 >;
    
    //! Default constructor
    NeighbourhoodManagerStrict()
    :particles{},centers{},particule_types{},unique_types{},type2idx{},positions{},pbc{},number_of_neighbours_stride{},neighbours{},neighbour_shifts{},number_of_neighbours{0},Rijs{this->implementation()},distances{this->implementation()}
    {}

    //! Copy constructor
    NeighbourhoodManagerStrict(const NeighbourhoodManagerStrict &other) = delete;

    //! Move constructor
    NeighbourhoodManagerStrict(NeighbourhoodManagerStrict &&other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerStrict() = default;

    //! Copy assignment operator
    NeighbourhoodManagerStrict& operator=(const NeighbourhoodManagerStrict &other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerStrict& operator=(NeighbourhoodManagerStrict &&other) = default;

    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    // return position vector
    inline Vector_ref get_position(const AtomRef_t& atom);

    // return number of center in the list
    inline size_t get_size() const;
    // return the id a given center atom
    inline size_t get_atom_id(const Parent& , int i_atom_id) const;

    //! return atom type
    inline int get_atom_type(const AtomRef_t& atom) const;

    // return the shift associated to the neighbour image. No implementation for Parent input type because centers are always within the cell
    template<int Level, int MaxLevel>
    inline Vector_ref get_atom_shift(const ClusterRef_t<Level, MaxLevel>& cluster, const int& j_cluster_id);


    // return the index of the center corresponding to its neighbour image
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster, int j_atom_id) const;

    // return the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

    void build(const Eigen::Ref<const Eigen::MatrixXd> positions, const Eigen::Ref<const VecXi>  particule_types,
               const Eigen::Ref<const VecXi> center_ids, 
               const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);

    inline AtomVectorProp_t::reference get_Rij(const ClusterRef_t<3,3>& cluster);

    inline AtomScalarProp_t::reference get_distance(const ClusterRef_t<3,3>& cluster);

    void update(const Eigen::Ref<const Eigen::MatrixXd> positions,const Eigen::Ref<const VecXi>  particule_types,
               const Eigen::Ref<const VecXi> center_ids, 
                const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);

  protected:

    void set_positions(const Eigen::Ref<const Eigen::MatrixXd> pos){
      this->positions = pos;
    }

    std::vector<AtomRef_t> particles;
    std::vector<AtomRef_t> centers; 
    std::vector<size_t> particule_types;
    std::vector<std::vector<AtomRef_t>> unique_types;
    std::map<int,int> type2idx;
    Matrix3XdC positions; 
    std::array<bool,3> pbc;
    std::vector<std::vector<size_t>> number_of_neighbours_stride;
    std::vector<std::vector<std::vector<AtomRef_t>>> neighbours;
    std::vector<std::vector<std::vector<Vector_t,Eigen::aligned_allocator<Vector_t>>>> neighbour_shifts;
    int number_of_neighbours;
    AtomVectorProp_t Rijs;
    AtomScalarProp_t distances;
    
  private:
  };

  /* ---------------------------------------------------------------------- */
  
  inline NeighbourhoodManagerStrict::AtomVectorProp_t::reference NeighbourhoodManagerStrict::get_Rij(
            const NeighbourhoodManagerStrict::ClusterRef_t<3,3>& cluster){
    return this->Rijs[cluster];
  }


  inline NeighbourhoodManagerStrict::AtomScalarProp_t::reference NeighbourhoodManagerStrict::get_distance(
            const NeighbourhoodManagerStrict::ClusterRef_t<3,3>& cluster){
    return this->distances[cluster];
  }
  /* ---------------------------------------------------------------------- */
  
  // return position vector
  inline Vector_ref NeighbourhoodManagerStrict::get_position(const NeighbourhoodManagerStrict::AtomRef_t& atom) {
    auto index{atom.get_index()};
    auto * xval{this->positions.col(index).data()};
    return Vector_ref(xval);
  }

  // return number of center in the list
  inline size_t NeighbourhoodManagerStrict::get_size() const {
    return this->centers.size();
  }
  // return the id a given center atom
  inline size_t NeighbourhoodManagerStrict::get_atom_id(const Parent& , int i_atom_id) const {
    return this->centers[i_atom_id].get_index();
  }

  //! return atom type
  inline int NeighbourhoodManagerStrict::get_atom_type(const NeighbourhoodManagerStrict::AtomRef_t& atom) const {
    auto && index{atom.get_index()};
    return this->particule_types[index];
  }
  
  template<int Level, int MaxLevel>
  inline Vector_ref NeighbourhoodManagerStrict::get_atom_shift(const NeighbourhoodManagerStrict::ClusterRef_t<Level, MaxLevel>& cluster, const int& j_cluster_id){
    auto && i_atom_id{cluster.get_atoms().front().get_index()};
    auto && type_idx{this->type2idx.at(cluster.get_atom_type())};
    auto * xval {this->neighbour_shifts[i_atom_id][type_idx][j_cluster_id].col(0).data()};
    return Vector_ref(xval);
  }


  // return the index of the center corresponding to its neighbour image
  template<>
  inline size_t NeighbourhoodManagerStrict::get_atom_id(const NeighbourhoodManagerStrict::ClusterRef_t<1,3>& cluster, int j_atom_id) const {
    auto && i_atom_id{cluster.get_atoms().front().get_index()};
    auto && ij_atom_id{this->unique_types[i_atom_id][j_atom_id].get_index()};
    return ij_atom_id;
  }

  template<>
  inline size_t NeighbourhoodManagerStrict::get_atom_id(const NeighbourhoodManagerStrict::ClusterRef_t<2,3>& cluster, int j_atom_id) const {
    auto && i_atom_id{cluster.get_atoms().front().get_index()};
    auto && type_idx{this->type2idx.at(cluster.get_atom_type())};
    auto && ij_atom_id{this->neighbours[i_atom_id][type_idx][j_atom_id].get_index()};
    return ij_atom_id;
  }
  
  // return the number of neighbours of a given atom
  template<>
  inline size_t NeighbourhoodManagerStrict::get_cluster_size(const NeighbourhoodManagerStrict::ClusterRef_t<1,3>& cluster) const {
    auto && i_atom_id{cluster.get_atoms().front().get_index()};
    size_t size{this->unique_types[i_atom_id].size()};
    return size;
  }

  // return the number of neighbours of a given atom
  template<>
  inline size_t NeighbourhoodManagerStrict::get_cluster_size(const NeighbourhoodManagerStrict::ClusterRef_t<2,3>& cluster) const {
    auto && i_atom_id{cluster.get_atoms().front().get_index()};
    auto && type_idx{this->type2idx.at(cluster.get_atom_type())};
    auto && size{this->neighbours[i_atom_id][type_idx].size()};
    return size;
  }
  
  template<>
  inline size_t NeighbourhoodManagerStrict::get_cluster_size(const NeighbourhoodManagerStrict::ClusterRef_t<3,3>& cluster) const {
    auto && i_atom_id{cluster.get_atoms().front().get_index()};
    auto && type_idx{this->type2idx.at(cluster.get_atom_type())};
    auto && size{this->neighbours[i_atom_id][type_idx].size()};
    return size;
  }

  /* ---------------------------------------------------------------------- */

  template<int Level, int MaxLevel>
  inline int NeighbourhoodManagerStrict::get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const {
    static_assert(Level == 3, "This class cas only handle single atoms and pairs");
    static_assert(MaxLevel == NeighbourhoodManagerStrict::traits::MaxLevel, "Wrong maxlevel");

    auto && icenter{cluster.get_atoms().front().get_index()};
    auto && type_idx{this->type2idx.at(cluster.get_atom_type())};
    auto && stride{this->number_of_neighbours_stride[icenter][type_idx]};
    auto && j{cluster.get_index()};
    auto main_offset{stride+j};
    return main_offset;
  }

  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms

  template <>
  inline int NeighbourhoodManagerStrict:: template get_offset_impl<1, 3>(const ClusterRef_t<1, 3>& cluster) const {
    return cluster.get_atoms().back().get_index();
  }
  /* ---------------------------------------------------------------------- */

    

}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_STRICT_H */
