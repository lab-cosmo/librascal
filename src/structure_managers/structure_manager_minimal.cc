/**
 * file   structure_manager_Minimal.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Implementation of Minimal neighbourhood manager
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

#include "structure_managers/structure_manager_minimal.hh"
#include "structure_managers/structure_manager_cell.hh"

namespace rascal {


  /* StructureManagerMinimal::Box StructureManagerMinimal::get_box(const int& bin_id){
    return this->boxes[bin_id];
  } */
  size_t StructureManagerMinimal::get_box_nb(){
    return this->boxes.size();
  }


   /* ---------------------------------------------------------------------- */

  void StructureManagerMinimal::build(const Eigen::Ref<const Eigen::MatrixXd>  positions,
                                        const Eigen::Ref<const VectorXi>  particule_types,
                                        const Eigen::Ref<const VectorXi> center_ids,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
    {
      Eigen::Index Natom{positions.cols()};

      const int dim{traits::Dim};

      // set the positions of all particle in the cell
      this->set_positions(positions);
      // set the references to the center positions
      for (int id{0}; id < center_ids.size(); ++id) {
          this->centers.push_back(StructureManagerMinimal::AtomRef_t(this->get_manager(),center_ids(id)));
      }
      // TODO get particles type as input and use it
      this->particule_types.resize(Natom);
      //set the references to the particles positions
      for (Eigen::Index id{0}; id < Natom; ++id){
        this->particles.push_back(StructureManagerMinimal::AtomRef_t(this->get_manager(),id));
        this->particule_types[id] = particule_types(id);
      }

      Cell_t lat = cell;
      this->lattice.set_cell(lat);
      this->cell = lat;
      Vec3_t reciprocal_lenghts = this->lattice.get_reciprocal_lenghts();
      double bin_size{cutoff_max};
      Vec3i_t nbins_c,neigh_search;
      Vec3_t nbins_cd;
      int nbins{1};
      double face_dist_c;


      for (int ii{0};ii<dim;++ii){
        // compute the distance between the cell faces (only French wiki https://fr.wikipedia.org/wiki/Distance_interr%C3%A9ticulaire)
        if (reciprocal_lenghts[ii] > 0){
          face_dist_c = 1 / reciprocal_lenghts[ii];
        }
        else {
          face_dist_c = 1;
        }
        // number of bin in each directions
        nbins_c[ii] =  std::max( static_cast<int>(face_dist_c/bin_size), 1);
        nbins_cd[ii] = static_cast<double>(nbins_c[ii]);
        // number of bin one need to look around
        neigh_search[ii] = static_cast<int>(std::ceil(bin_size * nbins_c[ii] / face_dist_c));
        // total number of bin
        nbins *= nbins_c[ii];

      }

      Vec3i_t bin_index_c;
      for (int ii{0}; ii < nbins; ++ii){
        internal::lin2mult<dim>(ii,nbins_c,bin_index_c);
        this->boxes.push_back(Box(this->get_manager(),bin_index_c, pbc, neigh_search, nbins_c));
      }

      // bin the particles in the boxes
      Vec3_t position_sc;
      int bin_id{0};
      this->part2bin.resize(Natom);
      for (auto part : this->particles){
          this->lattice.get_cartesian2scaled(part.get_position(),position_sc);
          bin_index_c = (position_sc.array() * nbins_cd.array()).cast<int>();
          bin_id = internal::mult2lin<dim>(bin_index_c,nbins_c);
          this->boxes[bin_id].push_particle_back(part.get_index());
          this->part2bin[part.get_index()] = bin_id;
      }

      // Set up the data strucure containing the information about neighbourhood
      // get the number of particles in the box and its neighbour boxes
      // set the arrays that will be used to iterate over the centers and neighbours
      this->neighbour_bin_id.resize(nbins);
      this->neighbour_atom_index.resize(nbins);
      //loop over the boxes
      for (size_t bin_index{0}; bin_index < this->boxes.size(); ++bin_index){
        size_t n_neigh{0};
        // loop over the neighbouring boxes
        for (size_t neigh_bin_id{0}; neigh_bin_id < this->boxes[bin_index].get_number_of_neighbour_box(); ++neigh_bin_id){
          int neig_bin_index{this->boxes[bin_index].get_neighbour_bin_index(neigh_bin_id)};
          //loop over the particle in the neighbouring boxes
          for (size_t neigh_part_id{0}; neigh_part_id < this->boxes[neig_bin_index].get_number_of_particles(); ++neigh_part_id){
            // store the indices to the corresponding atomic shift
            this->neighbour_bin_id[bin_index].push_back(AtomRef_t(this->get_manager(),neigh_bin_id));
            // store the indices to the neighbour particles
            this->neighbour_atom_index[bin_index].push_back(AtomRef_t(this->get_manager(),this->boxes[neig_bin_index].get_particle_index(neigh_part_id)));
          }
          n_neigh += this->boxes[neig_bin_index].get_number_of_particles();
        }
        this->boxes[bin_index].set_number_of_neighbours(n_neigh);
      }

      int stride{0};
      // get the stride for the fields. (center,neigh) dimentions are flattened in fields with center being leading dimension
      for (auto center : this->centers){
        this->number_of_neighbours_stride.push_back(stride);
        int bin_index{this->part2bin[center.get_index()]};
        size_t n_neigh{this->boxes[bin_index].get_number_of_neighbours()};
        stride += n_neigh;
      }
      this->number_of_neighbours = stride;
    }

  /* ---------------------------------------------------------------------- */

  StructureManagerMinimal::Box::Box(Manager_t& manager ,const Vec3i_t& coord, const std::array<bool, 3>& pbc,
        const Vec3i_t& neigh_search,const Vec3i_t& nbins_c)
        //const std::array<std::array<Dim_t, 3>,2>& neigh_bounds,
        :manager{manager},particles{},neighbour_bin_shift{},neighbour_bin_index{},number_of_neighbours{0},coordinates{}
  {
    const int dim{StructureManagerMinimal::dim()};

    this->coordinates = coord;
    int bin_id{0};
    Vec3i_t shift,neighbour_bin_idx_c;
    Vector_t neighbour_bin_shift;
    std::array<int,2> div_mod;
    std::array<std::vector<int>,dim> neigh_search_ids;

    // takes into account the pbc for neighbour boxes
    // TODO find a way to not have if statements
    for (int ii{0}; ii < dim; ++ii)  {
      if (pbc[ii] == true){
        for (int jj{-neigh_search[ii]}; jj <= neigh_search[ii]; ++jj ){
            neigh_search_ids[ii].push_back(this->coordinates[ii]+jj);
        }
      }
      else if (pbc[ii] == false && this->coordinates[ii] == 0 ){
        for (int jj{0}; jj <= neigh_search[ii]; ++jj ){
          neigh_search_ids[ii].push_back(this->coordinates[ii]+jj);
        }
      }
      else if (pbc[ii] == false && this->coordinates[ii] == nbins_c[ii]-1 ){
          for (int jj{-neigh_search[ii]}; jj <= 0; ++jj ){
          neigh_search_ids[ii].push_back(this->coordinates[ii]+jj);
        }
      }
    }

    for (auto dx : neigh_search_ids[0]){
      for (auto dy : neigh_search_ids[1]){
        for (auto dz : neigh_search_ids[2]){
          shift << dx,dy,dz;
          for (int ii{0};ii<dim;++ii){

            internal::div_mod(shift(ii),nbins_c(ii),div_mod);
            //internal::branchless_div_mod(shift(ii),nbins_c(ii),div_mod);
            neighbour_bin_shift[ii] = static_cast<double>(div_mod[0]);
            neighbour_bin_idx_c[ii] = div_mod[1];
          }

          bin_id = internal::mult2lin<dim>(neighbour_bin_idx_c,nbins_c);
          this->neighbour_bin_index.push_back(bin_id);
          this->neighbour_bin_shift.push_back(neighbour_bin_shift);
        }
      }
    }

  }

  inline size_t StructureManagerMinimal::Box::get_number_of_neighbour_box(){
    return this->neighbour_bin_index.size();
  }

  inline void StructureManagerMinimal::Box::push_particle_back(const int& part_index){
    this->particles.push_back(AtomRef_t(this->manager, part_index));
  }

  inline size_t StructureManagerMinimal::Box::get_number_of_particles(){
    return this->particles.size();
  }

  inline size_t StructureManagerMinimal::Box::get_number_of_neighbours(){
    return number_of_neighbours;
  }

  inline int StructureManagerMinimal::Box::get_neighbour_bin_index(const int& j_index){
    return this->neighbour_bin_index[j_index];
  }

  inline size_t StructureManagerMinimal::Box::get_particle_index(const int& index){
    return this->particles[index].get_index();
  }

  inline void StructureManagerMinimal::Box::set_number_of_neighbours(const size_t& neigh_nb){
    this->number_of_neighbours = neigh_nb;
  }


  /* ---------------------------------------------------------------------- */

  void StructureManagerMinimal::update(const Eigen::Ref<const Eigen::MatrixXd> positions,
                                        const Eigen::Ref<const VectorXi>  particule_types,
                                       const Eigen::Ref<const VectorXi> center_ids,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
  {
    bool some_condition{true};
    if (some_condition){
      StructureManagerMinimal::build(positions,particule_types,center_ids,cell,pbc,cutoff_max);
    }
  }

  /* ---------------------------------------------------------------------- */

  size_t StructureManagerMinimal::
  get_nb_clusters(int cluster_size) const {
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
