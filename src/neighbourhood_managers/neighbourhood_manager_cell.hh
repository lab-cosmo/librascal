/**
 * file   neighbourhood_manager_lammps.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager for lammps neighbourhood lists
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * proteus is distributed in the hope that it will be useful, but
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
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

namespace proteus {
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
    using Vector_block = typename Parent::Vector_block;
    using Vector_t = typename Parent::Vector_t;
    using AtomRef_t = typename Parent::AtomRef;
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template ClusterRef<Level, MaxLevel>;

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
    inline Vector_block get_position(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      //auto * xval{this->_positions.col(index)};
      //const Vector_t ppp = this->positions.col(index); //there is a copy I think
      //const Vector_ref aaa = Vector_ref(ppp.data());
      //Vector_ref aaa = Vector_ref(this->positions.col(index).data());
      //return aaa;
      return this->positions.col(index);
    }
    // return number of I atoms in the list
    inline Eigen::Index get_size() const {
      return this->positions.cols();
    }
    // return the number of neighbours of a given atom
    inline size_t get_atom_id(const Parent& , int i_atom_id) const {
      return this->centers[i_atom_id].get_index();
    }
    
    // return the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level, MaxLevel>& cluster) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      return this->neighlist[cluster.get_atoms().back().get_index()].size();
    }

    // return the number of atoms forming the next higher cluster with this one
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      return this->neighlist[std::move(i_atom_id)][j_atom_id].get_index();
    }


    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

    void build(const Eigen::MatrixXd& positions, const Eigen::MatrixXd& cell,const std::array<bool,3>& pbc, const double& cutoff_max);
    
    

  protected:
    std::vector<AtomRef_t> centers;
    std::vector<std::vector<AtomRef_t>> neighlist;
    std::vector<std::vector<std::array<int,3>>> offsetlist;
    Eigen::MatrixXd positions, cell;
    std::array<bool,3> pbc;
  private:
  };

  /* ---------------------------------------------------------------------- */
  void NeighbourhoodManagerCell::build(const Eigen::MatrixXd& positions,
                                        const Eigen::MatrixXd& cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
    {
      Eigen::Index Natom{positions.cols()};
      cout << "Natom " << Natom << endl;

      std::vector<AtomRef_t> particules;

      cout << "init centers " << endl;
      for (int id{0} ; id<Natom ; ++id) {
          this->centers.push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
          //cout << this->centers.at(id).get_position() << endl;
      }
      std::vector<std::vector<AtomRef_t>> neighlist(Natom);
      std::vector<std::vector<std::array<int,3>>> offsetlist(Natom);

      this->positions.resize(traits::Dim,Natom);
      this->cell.resize(traits::Dim,traits::Dim);

      this->positions = positions;
      this->cell = cell;
      this->pbc = pbc;

      std::vector<std::vector<int>> periodic_image_it(traits::Dim);
      for (unsigned int ii{0} ; ii < pbc.size() ; ++ii){
        if ( pbc[ii] ){
          int arr[3] = {-1,0,1};
          periodic_image_it[ii].assign(3,*arr);
        }
        else {
          int arr[1] = {0};
          periodic_image_it[ii].assign(1,*arr);
        }
      }

      cout << "pos " << endl << this->positions << endl;
      cout << "cell " << endl << this->cell << endl;
      //cout << "cell " << endl << this->cell << endl;
      Eigen::Array<double, traits::Dim, 1> cell_lengths;
      cell_lengths = cell.colwise().norm();
      cout << "cell_lengths " << endl << cell_lengths << endl;

      double bin_size{std::max(cutoff_max,3.)};
      cout << "bin_size " << endl << bin_size << endl;

      Eigen::Array<double, traits::Dim, 1> nbins_coord(0,0,0);
      nbins_coord = (cell_lengths / bin_size).ceil();
      cout << "nbins_coord " << endl << nbins_coord << endl;

      Eigen::Array<int, traits::Dim, 1> center_bin_coord(0,0,0);
      Eigen::Array<int, traits::Dim, 1> neigh_bin_coord(0,0,0);
      
      Eigen::Array<int, traits::Dim, 1> neigh_search(1,1,1);
      neigh_search = (bin_size * nbins_coord / cell_lengths).ceil().cast<int>();
      cout << "neigh_search " << endl << neigh_search << endl;

      Eigen::Array<int, traits::Dim, 1> upper_bound(1,1,1);
      Eigen::Array<int, traits::Dim, 1> lower_bound(1,1,1);
      double dist2{0};
      double cutoff_max2{cutoff_max*cutoff_max};
      Eigen::Vector3d offset(0,0,0);
      for (auto center: this->get_manager() ){
        // loop over the centers
        cout << "Index and positions " << center.get_index() << endl;
        center_bin_coord = (center.get_position().array() / nbins_coord).floor().cast<int>();
        cout << "center_bin_coord " << endl << center_bin_coord.transpose() << endl;
        upper_bound = center_bin_coord + neigh_search;
        lower_bound = center_bin_coord - neigh_search;
        std::vector<int> neighbours_index;
        std::vector<std::array<int,3> > offsets;
        for (int neigh_id{0} ; neigh_id<Natom ; ++neigh_id){ // TODO here particules is the same set of centers but this could be something else
        //for (auto neigh: this->get_manager() ){ 
          //int neigh_id{neigh.get_index()}; 
          // loop over the periodic images of neigh
          for (auto ii : periodic_image_it[0]){
            for (auto jj : periodic_image_it[1]){
              for (auto kk : periodic_image_it[2]){
                offset << ii,jj,kk;
                neigh_bin_coord = ( (positions.col(neigh_id) + (offset.transpose() * cell).transpose() ).array() / nbins_coord).floor().cast<int>();
                
                if ( (neigh_bin_coord < upper_bound).all() and (neigh_bin_coord > lower_bound).all()) {
                  //Vector_t distVec = ;
                  dist2 = ( (positions.col(neigh_id) + (offset.transpose() * cell).transpose() ) - center.get_position() ).squaredNorm();
                  if (dist2 <  cutoff_max2){
                    neighbours_index.push_back(neigh_id);
                    std::array<int, 3> arr = {ii,jj,kk};
                    offsets.push_back(arr);
                  }
                }
              }
            } 
          }
        }
        for (auto id: neighbours_index){
          neighlist[center.get_index()].push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
        }
        for (auto off : offsets){
          offsetlist[center.get_index()].push_back(off);
        }
        
      }

      for (auto center:centers){
      cout << "Center id: " << center << endl;
      cout << "Neighbour ids: " ;
      
        for (auto neigh : neighlist.at(center)){
            cout << neigh << ", "<< endl;
        }
      cout <<  endl;
      }

      this->neighlist =  std::move(neighlist);
      this->offsetlist =  std::move(offsetlist);

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
    //auto main_offset{this->offsets[i]};
    return j;
  }
  
  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms
  
  template <>
  inline int NeighbourhoodManagerCell:: template
  get_offset_impl<1, 2>(const ClusterRef_t<1, 2>& cluster) const {
    return cluster.get_atoms().back().get_index();
  }
  
}  // proteus

#endif /* NEIGHBOURHOOD_MANAGER_LAMMPS_H */
