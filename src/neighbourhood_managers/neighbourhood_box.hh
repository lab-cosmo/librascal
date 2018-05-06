/**
 * file   neighbourhood_boc.hh
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


#ifndef NEIGHBOURHOOD_BOX_H
#define NEIGHBOURHOOD_BOX_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/field.hh"
#include "lattice.hh"
#include "basic_types.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>


namespace rascal {

  template <class ManagerImplementation>
  class Box {
  public:
    using AtomRef_t = typename NeighbourhoodManagerBase<ManagerImplementation>::AtomRef;
    //! Default constructor
    Box() = default;

    //! constructor
    Box(NeighbourhoodManagerBase<ManagerImplementation> & manager, const Vec3i_t& coord, 
                const std::array<std::array<Dim_t, 3>,2>& neigh_bounds, const Vec3i_t& nbins_c)
                :manager{manager}
    { 
      //! warning: works only with negative a if |a| < b
      std::function<void (int,int,std::array<int,2>)> branchless_div_mod = [](int a, int b,std::array<int,2> d) {d = {(a+b)/b,(a+b)%b};};
      
      this->coordinates = coord;
      Vec3i_t shift,neighbour_bin_idx_c,neighbour_bin_shift;
      std::array<int,2> div_mod;
      int bin_id{0};
      for (int dx{neigh_bounds[0][0]}; dx <= neigh_bounds[1][0]; ++dx){
        for (int dy{neigh_bounds[0][1]}; dy <= neigh_bounds[1][1]; ++dy){
          for (int dz{neigh_bounds[0][2]}; dz <= neigh_bounds[1][2]; ++dz){
            shift << dx,dy,dz;
            
            for (int ii{0};ii<3;++ii){
              branchless_div_mod(coord(ii)+shift(ii),nbins_c(ii),div_mod);
              neighbour_bin_idx_c[ii] = div_mod[1];
              neighbour_bin_shift[ii] = div_mod[0];
            }

            bin_id = internal::mult2lin(neighbour_bin_idx_c,nbins_c);
            this->neighbour_bin_ids.push_back(bin_id);
            this->neighbour_bin_shift.push_back(neighbour_bin_shift);
            
          }
        }
      }
    };
    //! copy constructor
    Box(const Box & other) = default;
    //! assignment operator
    Box & operator=(const Box & other) = default;
    virtual ~Box() = default;

    inline void push_center_back(const int& id){
      this->centers.push_back(AtomRef_t(this->manager,id));
    }
    inline void push_neighbour_back(const int& id){
      this->neighbour_ids.push_back(AtomRef_t(this->manager,id));
    }
    
    inline size_t get_number_of_centers(){
      return this->centers.size();
    }

    inline size_t get_number_of_neighbour(){
      size_t size{this->neighbour_ids.size()};
      return size;
    }

    inline std::vector<AtomRef_t> get_neighbour_ids(){
      return this->neighbour_ids;
    }

    inline  int get_neighbour_index(int j_index){
      return this->neighbour_bin_ids[j_index];
    }

    inline  std::vector<AtomRef_t> get_centers(){
      return this->centers;
    }
    inline  std::vector<Dim_t> get_neighbour_bin_ids(){
      return this->neighbour_bin_ids;
    }

  

  protected:
    Vec3i_t coordinates;
    std::vector<AtomRef_t> centers;
    NeighbourhoodManagerBase<ManagerImplementation> & manager;
    std::vector<Vec3i_t,Eigen::aligned_allocator<Vec3i_t>> neighbour_bin_shift;
    std::vector<Dim_t> neighbour_bin_ids;
    // TODO replace neighbour_ids by an iterator that iterates over the centers in the neighbouring boxes
    std::vector<AtomRef_t> neighbour_ids; 
  };


  //----------------------------------------------------------------------------//
}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_CELL_H */
