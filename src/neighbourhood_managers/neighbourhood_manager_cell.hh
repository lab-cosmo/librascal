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
#include <izip.hh>

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
    using vVector3d = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >;

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
      //cout << "get_position" << endl;
      auto index{atom.get_index()};
      auto * xval{this->positions[index]};
      return Vector_ref(xval);
      //return this->positions[index];
    }
    
    // return number of I atoms in the list
    inline size_t get_size() const {
      return this->centers.size();
    }
    // return the number of neighbours of a given atom
    inline size_t get_atom_id(const Parent& , int i_atom_id) const {
      return this->centers[i_atom_id].get_index();
    }
    
    void set_positions(const Eigen::MatrixXd& pos){
      int dim{this->dim()};
      int Natom{pos.cols()};

      this->positions = new double*[Natom];
      for (int i{0}; i < Natom; ++i){
        this->positions[i] = new double[dim];
        for (int j{0}; j < dim; ++j){
          this->positions[i][j] = pos(j,i);
        }
      }
    }

    void set_positions(const std::vector<array<double,3>>& pos){
      int dim{this->dim()};
      int Natom{pos.size()};

      this->positions = new double*[Natom];
      for (int i{0}; i < Natom; ++i){
        this->positions[i] = new double[dim];
        for (int j{0}; j < dim; ++j){
          this->positions[i][j] = pos[i][j];
        }
      }
    }

    void set_cell(const Eigen::MatrixXd& cell){
      int dim{this->dim()};
      this->cell = new double*[dim];
      for (int i{0}; i < dim; ++i){
        this->cell[i] = new double[dim];
        for (int j{0}; j < dim; ++j){
          this->cell[i][j] = cell(j,i);
        }
      }
    }

    // return the number of atoms forming the next higher cluster with this one
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      auto && ij_atom_id{this->neighlist[i_atom_id][j_atom_id].get_index()};
      return ij_atom_id;
    }


    // return the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level, MaxLevel>& cluster) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      return this->neighlist[cluster.get_atoms().back().get_index()].size();
    }

    


    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

    void build(const Eigen::MatrixXd& positions, const Eigen::MatrixXd& cell,const std::array<bool,3>& pbc, const double& cutoff_max);
    
    void update(const Eigen::MatrixXd& positions, const Eigen::MatrixXd& cell,const std::array<bool,3>& pbc, const double& cutoff_max);
    

  protected:
    std::vector<AtomRef_t> centers; //! list of center indices. after build it points to neighpos
    std::vector<std::vector<AtomRef_t>> neighlist; //! list of list of neighbour indices. fisrt dimension is in the same order as centers. it points to neighpos
    std::vector<std::array<double,3>> neighpos; //! list of positions for the neighboor list. contains positions of the centers and then the positions of the neighboors.

    double **positions; //! list of pointers. after build it points to neighpos
    double **cell; //! list of pointers to the cell
    
    std::array<bool,3> pbc;
  private:
  };

  /* ---------------------------------------------------------------------- */
  
  void NeighbourhoodManagerCell::build(const Eigen::MatrixXd& positions,
                                        const Eigen::MatrixXd& cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
    {
      Eigen::Index Natom{positions.cols()};
      const int Ncenter{Natom};
      //this->centers.resize(Ncenter);
      this->neighlist.resize(Ncenter);

      cout << "Natom " << Natom << endl;
      const int dim{traits::Dim};
      using Vecd = typename Eigen::Array<double, dim, 1>;
      using Veci = typename Eigen::Array<int, dim, 1>;
      using Mati = typename Eigen::Array<int, dim, 2>;

      cout << "init centers " << endl;
      int count{0};
      for (int id{0} ; id<Ncenter ; ++id) {
          this->centers.push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
          array<double,3> iptpp = {positions(0,id),positions(1,id),positions(2,id)};
          this->neighpos.push_back(iptpp);
          count += 1;
          //cout << this->centers.at(id).get_position() << endl;
      }
      
      //std::vector<NeighbourhoodManagerCell::vVector3d> offsetlist(Natom);

      this->set_positions(positions);
      this->set_cell(cell);

      std::vector<std::vector<int>> periodic_image_it(dim);

      Vecd cell_lengths;
      cell_lengths = cell.colwise().norm();

      double bin_size{std::max(cutoff_max,3.)};

      //const Vecd nbins_coord(0,0,0);
      Vecd nbins_coord = (cell_lengths / bin_size).ceil();
      array<int,dim> nbins_c = {nbins_coord[0],nbins_coord[1],nbins_coord[2]};
      
      //Veci neigh_search(1,1,1);
      const Veci neigh_search = (bin_size * nbins_coord / cell_lengths).ceil().cast<int>();
      
      double dist2{0};
      double cutoff_max2{cutoff_max*cutoff_max};
      Eigen::Matrix<double, dim, 1> offset(0,0,0);
      Eigen::Matrix<double, dim, 1> jpos(0,0,0);
      Eigen::Matrix<double, dim, 1> jr(0,0,0);
      Veci center_bin_coord(0,0,0);
      Veci neigh_bin_coord(0,0,0);
      Veci bin_id(0,0,0);
      Mati bin_boundaries; //! 0->left 1->right 2->bottom 3->top
      bin_boundaries <<  0,nbins_c[0],
                         0,nbins_c[1],
                         0,nbins_c[2];

      vector<vector<vector<vector<int>>>> bin2icenter;
      vector<vector<vector<vector<array<int,dim>>>>> bin2neighbin,bin2neighbin_shift;
      bin2icenter.resize(nbins_coord(0));
      bin2neighbin.resize(nbins_coord(0));
      bin2neighbin_shift.resize(nbins_coord(0));
      for ( int ii{0} ; ii < nbins_coord[0] ; ++ii){
        bin2icenter[ii].resize(nbins_coord(1));
        bin2neighbin[ii].resize(nbins_coord(1));
        bin2neighbin_shift[ii].resize(nbins_coord(1));
        for ( int jj{0} ; jj < nbins_coord[1] ; ++jj){
          bin2icenter[ii][jj].resize(nbins_coord(2));
          bin2neighbin[ii][jj].resize(nbins_coord(2));
          bin2neighbin_shift[ii][jj].resize(nbins_coord(2));
          for ( int kk{0} ; kk < nbins_coord[2] ; ++kk){
            for ( int dx{-neigh_search[0]} ; dx <= neigh_search[0] ; ++dx){
              for ( int dy{-neigh_search[1]} ; dy <= neigh_search[1] ; ++dy){
                for ( int dz{-neigh_search[2]} ; dz <= neigh_search[2] ; ++dz){
                  bin_id << ii+dx,jj+dy,kk+dz;
                  array<int,dim> shift = {0,0,0};
                  array<int,dim> bin_id_c = {ii+dx,jj+dy,kk+dz};
                  array<int,dim> bin_id_new = {0,0,0};
                  if ( (bin_id >= bin_boundaries.col(0) ).all() and (bin_id < bin_boundaries.col(1) ).all()  ){
                    bin2neighbin[ii][jj][kk].push_back(bin_id_c);
                    bin2neighbin_shift[ii][jj][kk].push_back(shift);
                  }
                  else {
                    bool compatible_with_pbc = true;
                    for (int it{0};it<dim;++it){
                      shift[it] = (int)bin_id_c[it] / nbins_c[it];
                      //! remainder has the same sign as dividend 
                      bin_id_new[it] = bin_id_c[it] % nbins_c[it]; 
                      if  (bin_id_new[it] < 0){
                        bin_id_new[it] += nbins_c[it];
                        shift[it] -= 1;
                      }
                      if (pbc[it]==false and shift[it]!=0){
                        compatible_with_pbc = false;
                      }
                    }
                    if (compatible_with_pbc){
                      bin2neighbin[ii][jj][kk].push_back(bin_id_new);
                      bin2neighbin_shift[ii][jj][kk].push_back(shift);
                    }


                  }
              
                }
              
              }
            }
          }
        }
      }

      for (auto center: this->centers ){
        center_bin_coord = (center.get_position().array()/nbins_coord).floor().cast<int>();
      
        bin2icenter[center_bin_coord[0]][center_bin_coord[1]][center_bin_coord[2]].push_back(center.get_index());
      }
      
      for ( int ii{0} ; ii < nbins_coord[0] ; ++ii){
        for ( int jj{0} ; jj < nbins_coord[1] ; ++jj){
          for ( int kk{0} ; kk < nbins_coord[2] ; ++kk){
            for (auto icenter:bin2icenter[ii][jj][kk]){
              zipfor(jbin_id,jshift eachin bin2neighbin[ii][jj][kk],bin2neighbin_shift[ii][jj][kk]){
                //vector<int> ooo = bin2icenter[jbin_id[0]][jbin_id[1]][jbin_id[2]];
                for (auto jneigh:bin2icenter[jbin_id[0]][jbin_id[1]][jbin_id[2]]){
                  //Eigen::Map<Eigen::Matrix<double,1,dim>> dd(jshift,dim);
                  offset << jshift[0],jshift[1],jshift[2];
                  jpos = positions.col(jneigh) + (offset.transpose() * cell).transpose();
                  jr = jpos - positions.col(icenter);
                  dist2 = jr.squaredNorm();
                  if (cutoff_max2 > dist2){
                    this->neighlist[icenter].push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),count));
                    array<double,3> iptpp = {jpos[0],jpos[1],jpos[2]};
                    this->neighpos.push_back(iptpp);
                    count += 1;

                  }
                }
              }
            }
            
          }
        }
      }
      this->set_positions(this->neighpos);
      //Eigen::Array<double, dim, Natom> aa; = (positions/nbins_coord).ceil().cast<int>();

      /*
      for (unsigned int ii{0} ; ii < pbc.size() ; ++ii){
        if ( pbc[ii] ){
          std::vector<int> v = {-1,0,1};
          periodic_image_it[ii].insert( periodic_image_it[ii].end(), v.begin(), v.end() );

          //periodic_image_it[ii].push_back(v);
        }
        else {
          std::vector<int> v = {0};
          periodic_image_it[ii].insert( periodic_image_it[ii].end(), v.begin(), v.end() );

          //periodic_image_it[ii].push_back(v);
        }

      }
      for (int ii{0}; ii < 3; ++ii){
        cout << "p image " <<  ii  ;
        for (const auto& i: periodic_image_it[ii])
            cout << ' ' <<  i ;
        cout <<  endl;
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
      for (auto center: centers ){
        // loop over the centers
        cout << "Index and positions " << center.get_index() << endl;
        center_bin_coord = (center.get_position().array() / nbins_coord).floor().cast<int>();
        cout << "center_bin_coord "  << center_bin_coord.transpose() << endl;
        upper_bound = center_bin_coord + neigh_search;
        lower_bound = center_bin_coord - neigh_search;
        std::vector<int> neighbours_index;
        //std::vector<std::array<int,3> > offsets;
        NeighbourhoodManagerCell::vVector3d  offsets;
        for (int neigh_id{0} ; neigh_id<Natom ; ++neigh_id){ // TODO here particules is the same set of centers but this could be something else
        //for (auto neigh: this->get_manager() ){ 
          //int neigh_id{neigh.get_index()}; 
          // loop over the periodic images of neigh
          for (auto ii : periodic_image_it[0]){
            for (auto jj : periodic_image_it[1]){
              for (auto kk : periodic_image_it[2]){

                offset << ii,jj,kk;
                //cout << "offset : " << ii<<" "<< jj<<" "<< kk << endl;

                neigh_bin_coord = ( (positions.col(neigh_id) + (offset.transpose() * cell).transpose() ).array() / nbins_coord).floor().cast<int>();
                //cout << "Atom bin : " << neigh_bin_coord.transpose() << endl;
                if ( (neigh_bin_coord < upper_bound).all() and (neigh_bin_coord > lower_bound).all()) {
                  //Vector_t distVec = ;
                  dist2 = ( (positions.col(neigh_id) + (offset.transpose() * cell).transpose() ) - center.get_position() ).squaredNorm();
                  //cout << "Distance : " << dist2 << endl;
                  if (dist2 <  cutoff_max2){
                    //cout << "Id of neighbour: " << neigh_id << endl;
                    
                    cout << "Neighbour ids: " << neigh_id << endl;
                    cout << "neigh_bin_coord "  << neigh_bin_coord.transpose() << endl;
                    cout << "Offset: " << offset.transpose() << endl;
                    neighbours_index.push_back(neigh_id);
                    //std::array<int, 3> arr = {ii,jj,kk};
                    //offsets.push_back(arr);
                    offsets.push_back(offset);
                  }
                }
              }
            } 
          }
        }
        //cout << "# of neighbours: " << neighbours_index.size() << endl;
        for (auto id: neighbours_index){
          neighlist[center.get_index()].push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
        }
        for (auto off : offsets){
          offsetlist[center.get_index()].push_back(off);
        }
        cout  << endl;
      }

      
      
      this->neighlist =  std::move(neighlist);


      this->offsetlist =  std::move(offsetlist);

      
      for (auto center:this->centers){
        int c_idx{center.get_index()};
        cout << "Center id: " << c_idx << endl;
        cout << "Neighbour ids: " ;
          for (auto neigh : this->neighlist[c_idx]){
              cout << neigh.get_index() << ", ";
          }
          cout <<  endl;
          cout << "Neighbour shift: " <<  endl;
          for (auto off : this->offsetlist[c_idx]){
              cout << off[0] << ", "<<off[1] << ", "<<off[2] <<  " ; ";
          }
        cout <<  endl;
      }
      */
    }
  
  /* ---------------------------------------------------------------------- */
  
  void NeighbourhoodManagerCell::update(const Eigen::MatrixXd& positions,
                                        const Eigen::MatrixXd& cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
  {
    bool some_condition{false};
    if (some_condition){
      NeighbourhoodManagerCell::build(positions,cell,pbc,cutoff_max);
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
  
}  // proteus

#endif /* NEIGHBOURHOOD_MANAGER_LAMMPS_H */
