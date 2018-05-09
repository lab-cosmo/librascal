/**
 * file   neighbourhood_manager_cell.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Implementation of the neighbourhood manager for lammps
 *        neighbourhood lists
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

#include "neighbourhood_managers/neighbourhood_manager_cell.hh"

namespace rascal {
  

  NeighbourhoodManagerCell::Box NeighbourhoodManagerCell::get_box(const int& bin_id){
    return this->boxes[bin_id];
  }
  size_t NeighbourhoodManagerCell::get_box_nb(){
    return this->boxes.size();
  }
   /* ---------------------------------------------------------------------- */

  void NeighbourhoodManagerCell::build(const Eigen::Ref<const Eigen::MatrixXd>  positions,
                                        const std::vector<int>& center_ids,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
    {
      Eigen::Index Natom{positions.cols()};
      
      const int dim{traits::Dim};

      this->set_positions(positions);

      for (int id:center_ids) {
          this->centers.push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
      }
      for (Eigen::Index id{0}; id < Natom; ++id){
        this->particles.push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
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
      for (auto p : pbc){
        p;
      }
      std::array<std::array<Dim_t, 3>,2> neigh_bounds{{ {-neigh_search[0],-neigh_search[1],-neigh_search[2]},
                                                        { neigh_search[0], neigh_search[1], neigh_search[2]} }};
      for (int ii{0}; ii < nbins; ++ii){
        internal::lin2mult(ii,nbins_c,bin_index_c);
        this->boxes.push_back(Box(this->get_manager(),bin_index_c, neigh_bounds, nbins_c));
      }

      // bin the atoms in the boxes
      Vec3_t position_sc;
      int bin_id{0};
      for (auto part : this->particles){
          this->lattice.get_cartesian2scaled(part.get_position(),position_sc);
          bin_index_c = (position_sc.array() * nbins_cd.array()).cast<int>();
          bin_id = internal::mult2lin(bin_index_c,nbins_c);
          this->boxes[bin_id].push_particle_back(part.get_index());
          this->part2bin[part.get_index()] = bin_id;
      }
      
      // get the number of particles in the box and its neighbour boxes 
      // set the arrays that will be used to iterate over the centers and neighbours
      this->neighbour_bin_id.resize(nbins);
      this->neighbour_atom_index.resize(nbins);
      for (size_t bin_index{0}; bin_index < this->boxes.size(); ++bin_index){
        size_t n_neigh{0};
        for (size_t neigh_bin_id{0}; neigh_bin_id < this->boxes[bin_index].get_number_of_neighbour_box(); ++neigh_bin_id){
          
          int neig_bin_index{this->boxes[bin_index].get_neighbour_bin_index(neigh_bin_id)};
          for (size_t neigh_part_id{0}; neigh_part_id < this->boxes[neig_bin_index].get_number_of_particles(); ++neigh_part_id){
            // store the indices to the atomic shift inside the box
            this->neighbour_bin_id[bin_index].push_back(AtomRef_t(this->get_manager(),neigh_bin_id));
            // store the indices to the neighbour particles in positions
            this->neighbour_atom_index[bin_index].push_back(AtomRef_t(this->get_manager(),this->boxes[neig_bin_index].get_particle_index(neigh_part_id)));
          }
          n_neigh += this->boxes[neig_bin_index].get_number_of_particles();
        }
        this->boxes[bin_index].set_number_of_neighbours(n_neigh);
      }

      int stride{0};
      for (auto center : this->centers){
        this->number_of_neighbours_stride.push_back(stride);
        int bin_index{this->part2bin[center.get_index()]};
        size_t n_neigh{this->boxes[bin_index].get_number_of_neighbours()};
        stride += n_neigh;
      }
      this->number_of_neighbours = stride;
    }

  /* ---------------------------------------------------------------------- */

  NeighbourhoodManagerCell::Box::Box(Manager_t& manager ,const Vec3i_t& coord,
        const std::array<std::array<Dim_t, 3>,2>& neigh_bounds, 
        const Vec3i_t& nbins_c)
        :manager{manager},particles{},neighbour_bin_shift{},neighbour_bin_index{},number_of_neighbours{0},coordinates{}
  { 
    this->coordinates = coord;
    Vec3i_t shift,neighbour_bin_idx_c;
    Vector_t neighbour_bin_shift;
    std::array<int,2> div_mod;
    int bin_id{0};
    for (int dx{neigh_bounds[0][0]}; dx <= neigh_bounds[1][0]; ++dx){
      for (int dy{neigh_bounds[0][1]}; dy <= neigh_bounds[1][1]; ++dy){
        for (int dz{neigh_bounds[0][2]}; dz <= neigh_bounds[1][2]; ++dz){
          shift << dx,dy,dz;
          
          for (int ii{0};ii<3;++ii){
            //branchless_div_mod(coord(ii)+shift(ii),nbins_c(ii),div_mod);
            internal::div_mod(coord(ii)+shift(ii),nbins_c(ii),div_mod);
            neighbour_bin_shift[ii] = static_cast<double>(div_mod[0]);
            neighbour_bin_idx_c[ii] = div_mod[1];
          }

          bin_id = internal::mult2lin(neighbour_bin_idx_c,nbins_c);
          this->neighbour_bin_index.push_back(bin_id);
          this->neighbour_bin_shift.push_back(neighbour_bin_shift);
          
        }
      }
    }
  }

  inline size_t NeighbourhoodManagerCell::Box::get_number_of_neighbour_box(){
    return this->neighbour_bin_index.size();
  }

  inline void NeighbourhoodManagerCell::Box::push_particle_back(const int& part_index){
    this->particles.push_back(AtomRef_t(this->manager, part_index));
  }

  inline size_t NeighbourhoodManagerCell::Box::get_number_of_particles(){
    return this->particles.size();
  }
  
  inline size_t NeighbourhoodManagerCell::Box::get_number_of_neighbours(){
    return number_of_neighbours;
  }

  inline int NeighbourhoodManagerCell::Box::get_neighbour_bin_index(const int& j_index){
    return this->neighbour_bin_index[j_index];
  }

  inline size_t NeighbourhoodManagerCell::Box::get_particle_index(const int& index){
    return this->particles[index].get_index();
  }
  
  inline void NeighbourhoodManagerCell::Box::set_number_of_neighbours(const size_t& neigh_nb){
    this->number_of_neighbours = neigh_nb;
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

  size_t NeighbourhoodManagerCell::
  get_nb_clusters(int cluster_size)  {
    switch (cluster_size) {
    case 1: {
      return this->centers.size();
      break;
    }
    case 2: {
      return this->number_of_neighbours;
      break;
    }
    default:
      throw std::runtime_error("Can only handle single atoms and pairs");
      break;
    }
  }

}  // rascal
