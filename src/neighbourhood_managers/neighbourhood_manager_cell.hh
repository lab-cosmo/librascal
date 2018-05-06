/**
 * file   neighbourhood_manager_cell.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager for lammps neighbourhood lists
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


#ifndef NEIGHBOURHOOD_MANAGER_CELL_H
#define NEIGHBOURHOOD_MANAGER_CELL_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/neighbourhood_box.hh"
#include "neighbourhood_managers/field.hh"
#include "lattice.hh"
#include "basic_types.h"
#include "neighbourhood_managers/field.hh"
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
  class NeighbourhoodManagerCell;

  //! traits specialisation for Lammps manager
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerCell> {
    constexpr static int Dim{3};
    constexpr static int MaxLevel{2};
  };
  class NeighbourhoodManagerCell: public NeighbourhoodManagerBase<NeighbourhoodManagerCell>
  {
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerCell>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerCell>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using AtomRef_t = typename Parent::AtomRef;
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template ClusterRef<Level, MaxLevel>;
    using vVector3d = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >;
    using AtomVectorField_t = Field<NeighbourhoodManagerCell, double, 1, 3>;
    //! Default constructor
    NeighbourhoodManagerCell() = default;

    //! Copy constructor
    NeighbourhoodManagerCell(const NeighbourhoodManagerCell &other) = delete;

    //! Move constructor
    NeighbourhoodManagerCell(NeighbourhoodManagerCell &&other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerCell() = default;

    //! Copy assignment operator
    NeighbourhoodManagerCell& operator=(const NeighbourhoodManagerCell &other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerCell& operator=(NeighbourhoodManagerCell &&other) = default;

    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    // return position vector
    inline Vector_ref get_position(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto * xval{this->positions.col(index).data()};
      return Vector_ref(xval);
    }

    // return number of I atoms in the list
    inline size_t get_size() const {
      return this->centers.size();
    }
    // return the number of neighbours of a given atom
    inline size_t get_atom_id(const Parent& , int i_atom_id) const {
      return this->centers[i_atom_id].get_index();
    }

    void set_positions(const Eigen::Ref<const Eigen::MatrixXd> pos){
      this->positions = pos;
    }

    // return the number of atoms forming the next higher cluster with this one
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster, int j_atom_id) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      auto && i_bin_id{this->center2bin.at(i_atom_id)};
      Box<NeighbourhoodManagerCell> box = this->boxes[i_bin_id];
      auto && ij_atom_id{box.get_neighbour_index(j_atom_id)};
      return ij_atom_id;
    }


    // return the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level, MaxLevel>& cluster) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      auto && box_id{this->center2bin.at(cluster.get_atoms().back().get_index())};
      Box<NeighbourhoodManagerCell> box = this->boxes[box_id];
      int aa{box.get_number_of_neighbour()};
      return aa;
    }

    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

    size_t get_nb_of_center_in_box(const int& bin_id){
      return this->boxes[bin_id].get_number_of_centers();
    }

    std::vector<AtomRef_t> get_centers_from_box(int bin_id){
      Box<NeighbourhoodManagerCell> box = this->boxes[bin_id];
      return box.get_centers();
    }
    //void build(const Eigen::MatrixXd& positions, const Eigen::MatrixXd& cell,const std::array<bool,3>& pbc, const double& cutoff_max);
    void build(const Eigen::Ref<const Eigen::MatrixXd> positions,const std::vector<int>& center_ids, const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);


    void update(const Eigen::Ref<const Eigen::MatrixXd>,const std::vector<int>& center_ids, const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);


  protected:
    std::vector<AtomRef_t> centers; //! list of center indices. after build it points to neighpos
    std::vector<std::vector<AtomRef_t>> neighlist; //! list of list of neighbour indices. fisrt dimension is in the same order as centers. it points to neighpos
    std::vector<std::array<double,3>> neighpos; //! list of positions for the neighboor list. contains positions of the centers and then the positions of the neighboors.
    Matrix3XdC positions; //! list of pointers. after build it points to neighpos
    Lattice lattice;
    std::array<bool,3> pbc;
    std::vector<Box<NeighbourhoodManagerCell>> boxes;
    std::map<int, int> center2bin; 
  private:
  };

  /* ---------------------------------------------------------------------- */

  void NeighbourhoodManagerCell::build(const Eigen::Ref<const Eigen::MatrixXd>  positions,
                                        const std::vector<int>& center_ids,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
    {
      Eigen::Index Natom{positions.cols()};
      const int Ncenter{Natom};
      //this->centers.resize(Ncenter);
      this->neighlist.resize(Ncenter);

      const int dim{traits::Dim};

      int count{0};
      this->set_positions(positions);

      for (int id:center_ids) {
          this->centers.push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
          count += 1;
      }

      Cell_t lat = cell;
      this->lattice.set_cell(lat);
      Vec3_t reciprocal_lenghts = this->lattice.get_reciprocal_lenghts();
      double bin_size{cutoff_max};
      Vec3i_t nbins_c,neigh_search;
      Vec3_t nbins_cd;
      int nbins{1};
      double face_dist_c;
      for (int ii{0};ii<dim;++ii){
        if (reciprocal_lenghts[ii] > 0){
          face_dist_c = 1 / reciprocal_lenghts[ii];
        }
        else {
          face_dist_c = 1;
        }
      
        nbins_c[ii] =  std::max( static_cast<int>(face_dist_c/bin_size), 1);
        nbins_cd[ii] = static_cast<double>(nbins_c[ii]);
        neigh_search[ii] = static_cast<int>(std::ceil(bin_size * nbins_c[ii] / face_dist_c));
        nbins *= nbins_c[ii];
      
      }

      Vec3i_t bin_index_c;
      // TODO take into acount pbc dif from 1,1,1 
      std::array<std::array<Dim_t, 3>,2> neigh_bounds{{{-nbins_c[0],-nbins_c[1],-nbins_c[2]},{nbins_c[0],nbins_c[1],nbins_c[2]}}};
      for (int ii{0}; ii < nbins; ++ii){
        internal::lin2mult(ii,nbins_c,bin_index_c);
        this->boxes.push_back(Box<NeighbourhoodManagerCell>(this->get_manager(), bin_index_c, neigh_bounds, nbins_c));
      }

      Vec3_t position_sc;
      int bin_id{0};
      for (auto center : this->centers){
          this->lattice.get_cartesian2scaled(center.get_position(),position_sc);
          bin_index_c = (position_sc.array() * nbins_cd.array()).cast<int>();
          bin_id = internal::mult2lin(bin_index_c,nbins_c);
          this->boxes[bin_id].push_center_back(center.get_index());
          this->center2bin[center.get_index()] = bin_id;
      }

      for (auto box : this->boxes){
        for (auto neigh_bin : box.get_neighbour_bin_ids()) {
          for (auto neigh : this->boxes[neigh_bin].get_centers()){
            box.push_neighbour_back(neigh.get_index());
            size_t aa{box.get_number_of_neighbour()};
          }
        }
      }
      
    }

  /* ---------------------------------------------------------------------- */

  void NeighbourhoodManagerCell::update(const Eigen::Ref<const Eigen::MatrixXd> positions,
                                        const std::vector<int>& center_ids,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
  {
    bool some_condition{false};
    if (some_condition){
      NeighbourhoodManagerCell::build(positions,center_ids,cell,pbc,cutoff_max);
    }
  }
  /* ---------------------------------------------------------------------- */

  template<int Level, int MaxLevel>
  inline int NeighbourhoodManagerCell::
  get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const {
    static_assert(Level == 2, "This class cas only handle single atoms and pairs");
    static_assert(MaxLevel == traits::MaxLevel, "Wrong maxlevel");

    auto atoms{cluster.get_atoms()};
    auto i{atoms.front().get_index()};
    auto j{cluster.get_index()};
    auto main_offset{this->neighlist[i][j]};
    return main_offset;
  }

  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms

  template <>
  inline int NeighbourhoodManagerCell:: template
  get_offset_impl<1, 2>(const ClusterRef_t<1, 2>& cluster) const {
    return cluster.get_atoms().back().get_index();
  }
  /* ---------------------------------------------------------------------- */
  
  

  //----------------------------------------------------------------------------//
}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_CELL_H */
