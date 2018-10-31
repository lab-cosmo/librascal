/**
 * file   representation_manager_sorted_coulomb.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  base class for representation managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef BASIS_REPRESENTATION_MANAGER_SORTED_COULOMB_H
#define BASIS_REPRESENTATION_MANAGER_SORTED_COULOMB_H


#include "representations/representation_manager_base.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"
#include <vector>
#include <algorithm>
#include <string>


namespace rascal {

  namespace internal {

    typedef std::vector<double>::const_iterator myiter;

    struct ordering {
      bool operator ()(std::pair<size_t, myiter> const& a, 
                              std::pair<size_t, myiter> const& b) {
          return *(a.second) < *(b.second);
      }
    };


    template <typename T>
    std::vector<T> sort_from_ref(
        std::vector<T> const& in,
        std::vector<std::pair<size_t, myiter> > const& reference) {
        std::vector<T> ret(in.size());

        size_t const size = in.size();
        for (size_t i{0}; i < size; ++i)
          {
            ret[i] = in[reference[i].first];
          }

        return ret;
    }

    template <typename DerivedA,typename DerivedB>
    void sort_coulomb_matrix(
      const Eigen::DenseBase<DerivedA> & in,
      Eigen::DenseBase<DerivedB> & out,
      std::vector<std::pair<size_t, myiter> > const& reference) {
      
      size_t ncol = in.cols();
      size_t ii{0};
      // auto nn{ncol*(ncol-1)/2};
      // std::cout <<ncol<<", "<<nn <<std::endl;
      for (size_t i{0}; i < ncol; ++i) {
        auto col_id{reference[i].first};
        for (size_t j{col_id}; j < ncol; ++j) {
          out(ii) = in(j,col_id);
          ii += 1;
        }
      }
    }
  }

  template<class StructureManager>
  class RepresentationManagerSortedCoulomb: public RepresentationManagerBase
  {
  public:
    // TODO make a traits mechanism
    using hypers_t = RepresentationManagerBase::hypers_t;
    using Property_t = Property<double, 1, 1, Eigen::Dynamic, 1>;
    using Manager_t = StructureManager;

    //! Default constructor 
    RepresentationManagerSortedCoulomb(Manager_t &sm, 
      double central_decay , double interaction_cutoff, 
      double interaction_decay, size_t size)
      :structure_manager{sm},central_decay{central_decay},
      interaction_cutoff{interaction_cutoff},
      interaction_decay{interaction_decay},size{size},
      coulomb_matrices{sm}
      {
        for (auto center: this->structure_manager){
          auto Nneighbours{center.size()};
          if (Nneighbours > this->size){
              std::cout << "size is too small for this "
                           "structure and has been reset to: " 
                        << Nneighbours << std::endl;
              this->size = Nneighbours;
          }
        }
      }

    RepresentationManagerSortedCoulomb(Manager_t &sm,const hypers_t& hyper)
      :structure_manager{sm},central_decay{},
      interaction_cutoff{},
      interaction_decay{},coulomb_matrices{sm}
      {
        this->set_hyperparameters(hyper);
      }

    //! Copy constructor
    RepresentationManagerSortedCoulomb(
      const RepresentationManagerSortedCoulomb &other) = delete;

    //! Move constructor
    RepresentationManagerSortedCoulomb(
      RepresentationManagerSortedCoulomb &&other) = default;

    //! Destructor
    virtual ~RepresentationManagerSortedCoulomb()  = default;

    //! Copy assignment operator
    RepresentationManagerSortedCoulomb& operator=(
      const RepresentationManagerSortedCoulomb &other) = delete;

    //! Move assignment operator
    RepresentationManagerSortedCoulomb& operator=(
      RepresentationManagerSortedCoulomb && other) = default;
    

    //! compute representation
    void compute();

    //! getter for the representation
    Eigen::Map<Eigen::MatrixXd> get_representation_full() {
      auto Nb_centers{this->structure_manager.nb_clusters(1)};
      auto Nb_features{this->size};
      auto& raw_data{this->coulomb_matrices.get_raw_data()};
      Eigen::Map<Eigen::MatrixXd> representation(raw_data.data(),Nb_features,Nb_centers);
      std::cout <<"Sizes: "<< representation.size()<<", "<< Nb_features<<", "<<Nb_centers<<std::endl;
      return representation;
    }

    // TODO think of a generic input type for the hypers
    void set_hyperparameters(const hypers_t & );

    Manager_t& structure_manager;
    //hypers_t hyperparmeters;
    double central_decay;
    double interaction_cutoff;
    double interaction_decay;
    // first dimension of the largest coulomb mat
    size_t size;

    const std::string name = internal::GetTypeName<StructureManager>();

    Property_t coulomb_matrices;

  protected:
  private:
  };


  template<class Mngr>
  void RepresentationManagerSortedCoulomb<Mngr>::compute(){
    // initialise the coulomb_matrices storage

    this->coulomb_matrices.resize_to_zero();
    this->coulomb_matrices.set_nb_row(this->size);
    std::cout << "name is: " << this->name << std::endl;
    typedef std::vector<double>::const_iterator distiter;
    std::cout << this->structure_manager.nb_clusters(1)<< ", "<<
                  this->structure_manager.nb_clusters(2)  <<", "<<
                 this->structure_manager.size() <<std::endl;

    for (auto center: this->structure_manager){

      // auto pos = center.get_position();
      // for (auto jj{0}; jj < pos.size(); ++jj){
      //       std::cout << pos(jj) << ",\t";
      // }
      // std::cout << std::endl;
      std::vector<double> distances_to_sort{};
      // the local coulomb matrix
      auto Nneighbours{center.size()};
      Eigen::MatrixXd coulomb_mat(Nneighbours,Nneighbours);
      Eigen::MatrixXd lin_sorted_coulomb_mat(this->size*(this->size+1)/2,1);

      for (auto neigh_i: center){
        size_t ii{neigh_i.get_index()};
        auto dik{this->structure_manager.get_distance(neigh_i)};
        distances_to_sort.push_back(dik);
        auto Zi{neigh_i.get_atom_type()};
        coulomb_mat(ii,ii) = 0.5*std::pow(Zi,2.4);
        for (auto neigh_j: center){
          size_t jj{neigh_j.get_index()};
          // work only on the lower diagonal
          if (ii >= jj) continue;
          auto Zj{neigh_i.get_atom_type()};
          auto dij{(neigh_i.get_position()-neigh_j.get_position()).norm()};
          coulomb_mat(jj,ii) = Zi*Zj/dij;
        }
      }

      std::vector<std::pair<size_t, distiter> > order_coulomb(distances_to_sort.size());

      size_t nn{0};
      for (distiter it{distances_to_sort.begin()}; 
                            it != distances_to_sort.end(); ++it, ++nn)
          {order_coulomb[nn] = make_pair(nn, it);}
      
      std::sort(order_coulomb.begin(), order_coulomb.end(), internal::ordering());

      internal::sort_coulomb_matrix(coulomb_mat,lin_sorted_coulomb_mat,order_coulomb);
      this->coulomb_matrices.push_back(lin_sorted_coulomb_mat);
      // std::cout << "Compute atomic CM: "<< center.get_atom_index() << std::endl;
      // for (auto ii{0}; ii < coulomb_mat.cols(); ++ii){
      //   for (auto jj{0}; jj < coulomb_mat.rows(); ++jj){
      //       std::cout << coulomb_mat(jj,ii) << ",\t";
      //   }
      //   std::cout << std::endl;
      // }
      // std::cout << "Compute atomic CM: "<< center.get_atom_index() << std::endl;
      // for (auto jj{0}; jj < lin_sorted_coulomb_mat.rows(); ++jj){
      //           std::cout << lin_sorted_coulomb_mat(jj,0) << ", ";
      //       }
      //       std::cout << std::endl;
    }



  }

}

#endif /* BASIS_REPRESENTATION_MANAGER_SORTED_COULOMB_H */

