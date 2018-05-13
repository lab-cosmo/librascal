/**
 * file   neighbourhood_manager_strict.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Implementation of the linked cell neighbourhood manager
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

#include "neighbourhood_managers/neighbourhood_manager_strict.hh"

namespace rascal {

  /* ---------------------------------------------------------------------- */

  void NeighbourhoodManagerStrict::build(const Eigen::Ref<const Eigen::MatrixXd>  positions,
                                        const Eigen::Ref<const VecXi>  particule_types,
                                        const Eigen::Ref<const VecXi> center_ids,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
  { 
    using Vector_t = NeighbourhoodManagerStrict::Vector_t;
    using AtomRef_t = NeighbourhoodManagerStrict::AtomRef_t;
    Eigen::Index Natom{positions.cols()};
    
    //const int dim{traits::Dim};

    // set the positions of all particle in the cell
    this->set_positions(positions);
    Cell_t lattice = cell;
    // set the references to the center positions
    for (int id{0}; id < center_ids.size(); ++id) {
        this->centers.push_back(AtomRef_t(this->get_manager(),center_ids(id)));
    }
    // TODO get particles type as input and use it
    this->particule_types.resize(Natom);
    //set the references to the particles positions
    std::vector<int> unique_types;
    std::vector<int> unique_type_part_ids;
    for (Eigen::Index id{0}; id < Natom; ++id){
      this->particles.push_back(AtomRef_t(this->get_manager(),id));
      this->particule_types[id] = particule_types(id);
      if(std::find(unique_types.begin(), unique_types.end(), particule_types(id)) == unique_types.end()) {
          unique_types.push_back(particule_types(id));
          unique_type_part_ids.push_back(id);
      } 
      //unique_types.insert(particule_types(id));
    }
    //this->unique_types.assign( unique_types.begin(), unique_types.end() ); // unique_types is sorted

    for (size_t ii{0}; ii < unique_types.size(); ++ii){
      this->type2idx[unique_types[ii]] = ii;
      //this->unique_types.push_back(AtomRef_t(this->get_manager(),unique_type_part_ids[ii]));
    }

    NeighbourhoodManagerCell cell_manager;
    cell_manager.build(positions,particule_types,center_ids,cell,pbc,cutoff_max);

    size_t Ncenter{this->centers.size()};
    number_of_neighbours_stride.resize(Ncenter);
    this->neighbours.resize(Ncenter);
    this->unique_types.resize(Ncenter);
    this->neighbour_shifts.resize(Ncenter);
    this->distances.resize_to_zero();
    this->Rijs.resize_to_zero();

    double rc2{cutoff_max*cutoff_max};
    this->number_of_neighbours = 0;
    for (auto center: cell_manager) {
      int center_idx{center.get_atom_index()};
      size_t nb_types{unique_types.size()};
      number_of_neighbours_stride[center_idx].resize(nb_types);
      this->neighbours[center_idx].resize(nb_types);
      this->neighbour_shifts[center_idx].resize(nb_types);
      std::vector<std::vector<Vector_t,Eigen::aligned_allocator<Vector_t>> > Rijs(nb_types);
      std::vector<std::vector<double>> distances(nb_types);
      std::vector<size_t> nb_of_neigh_inner(nb_types,0);
      int n_neigh_outer{this->number_of_neighbours};

      for (auto neigh: center) {
        double d2{(neigh.get_position() + lattice * neigh.get_atom_shift() - center.get_position()).squaredNorm()};
        if (d2 < rc2){
          int type_idx{this->type2idx.at(neigh.get_atom_type())};


          Vector_t Rij = neigh.get_position() + lattice * neigh.get_atom_shift() - center.get_position();
          this->neighbours[center_idx][type_idx].push_back(AtomRef_t(this->get_manager(),neigh.get_atom_index()));
          Vector_t shift = neigh.get_atom_shift();
          this->neighbour_shifts[center_idx][type_idx].push_back(shift);
          Rijs[type_idx].push_back(Rij);
          distances[type_idx].push_back(std::sqrt(d2));
          nb_of_neigh_inner[type_idx] += 1;
          ++this->number_of_neighbours;
        }      
      }

      int stride{0};
      for (size_t type_idx{0}; type_idx < nb_types; ++type_idx){

        
        for (size_t jneigh{0}; jneigh < this->neighbours[center_idx][type_idx].size(); ++jneigh){
          this->distances.push_back(distances[type_idx][jneigh]);
          //auto * xval{Rijs[type_idx][jneigh].col(0).data()} ; 
          Eigen::Map<Vector_t> v(&Rijs[type_idx][jneigh](0),Rijs[type_idx][jneigh].size());
          this->Rijs.push_back(v);
        }

        if (this->neighbours[center_idx][type_idx].size() > 0){
          this->unique_types[center_idx].push_back(AtomRef_t(this->get_manager(),unique_type_part_ids[type_idx]));
          this->number_of_neighbours_stride[center_idx][type_idx] = n_neigh_outer + stride;
          stride += nb_of_neigh_inner[type_idx];
        }
        
      }
      
    }
 

  }


  /* ---------------------------------------------------------------------- */

  void NeighbourhoodManagerStrict::update(const Eigen::Ref<const Eigen::MatrixXd> positions,
                                        const Eigen::Ref<const VecXi>  particule_types,
                                        const Eigen::Ref<const VecXi> center_ids,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
  {
    bool some_condition{false};
    if (some_condition){
      NeighbourhoodManagerStrict::build(positions,particule_types,center_ids,cell,pbc,cutoff_max);
    }
  }

  /* ---------------------------------------------------------------------- */

  size_t NeighbourhoodManagerStrict::get_nb_clusters(int cluster_size)  {
      switch (cluster_size) {
      case 1: {
        return this->centers.size();
        break;
      }
      /* cluster_size == 2 is a dummy dimension ordering the neighbour iteration by species
      case 2: {
        return this->unique_types.size();
        break;
      }*/
      case 3: {
        return this->number_of_neighbours;
        break;
      }
      default:
        throw std::runtime_error("Can only handle single atoms and pairs");
        break;
      }
    }

}  // rascal
