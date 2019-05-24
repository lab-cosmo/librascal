/**
 * file   representation_manager_sorted_coulomb.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  Implements the Sorted Coulomb representation
 *
 * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SORTED_COULOMB_HH_
#define SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SORTED_COULOMB_HH_

#include "representations/representation_manager_base.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"

#include <vector>
#include <algorithm>
#include <math.h>

namespace rascal {

  namespace internal {

    using distiter = typename std::vector<double>::const_iterator;

    /**
     * Function for the sorting of a container using the order
     * from another container
     */
    struct ordering {
      static bool ascending(std::pair<size_t, distiter> const & a,
                            std::pair<size_t, distiter> const & b) {
        return *(a.second) < *(b.second);
      }

      static bool descending(std::pair<size_t, distiter> const & a,
                             std::pair<size_t, distiter> const & b) {
        return *(a.second) > *(b.second);
      }
    };

    /* ---------------------------------------------------------------------- */

    /* -------------------- rep-options-def-start -------------------- */
    // Enum class defining the several possible sorting options of the Coulomb
    // Matrix
    enum class CMSortAlgorithm {
      Distance,  // sort according to the distance from the central atom
      RowNorm,   // sort according to the norm of the coulomb matrix rows
    };

    // Empty general template class implementing the determination of the
    // sorting order for the coulomb matrix. It should never be instantiated.
    template <CMSortAlgorithm Method>
    struct SortCoulomMatrix {};
    /* -------------------- rep-options-def-end -------------------- */

    /* -------------------- rep-options-impl-start -------------------- */
    /**
     * Sort the coulomb matrix using the distance to the central atom
     * as reference order.
     *
     * @params distance_mat distance matrix between all the atoms in the
     *                      neighbourhood
     */
    template <>
    struct SortCoulomMatrix<CMSortAlgorithm::Distance> {
      static decltype(auto) get_coulomb_matrix_sorting_order(
          const Eigen::Ref<const Eigen::MatrixXd> & distance_mat,
          const Eigen::Ref<const Eigen::MatrixXd> &) {
        // initialize the distances to be sorted. the center is always first
        std::vector<double> distances_to_sort{0};
        distances_to_sort.reserve(distance_mat.cols());

        for (auto idx_i{1}; idx_i < distance_mat.cols(); ++idx_i) {
          distances_to_sort.push_back(distance_mat(idx_i, 0));
        }

        // find the sorting order
        std::vector<std::pair<size_t, distiter>> order_coulomb(
            distances_to_sort.size());
        size_t nn{0};
        for (distiter it{distances_to_sort.begin()};
             it != distances_to_sort.end(); ++it, ++nn) {
          order_coulomb[nn] = make_pair(nn, it);
        }

        // use stable sort
        std::stable_sort(order_coulomb.begin(), order_coulomb.end(),
                         ordering::ascending);
        return order_coulomb;
      }
    };

    /**
     * Sort the coulomb matrix using the distance to the central atom
     * as reference order.
     *
     * @params coulomb_mat coulomb matris between all the atoms in the
     *                      neighbourhood
     */
    template <>
    struct SortCoulomMatrix<CMSortAlgorithm::RowNorm> {
      static decltype(auto) get_coulomb_matrix_sorting_order(
          const Eigen::Ref<const Eigen::MatrixXd> &,
          const Eigen::Ref<const Eigen::MatrixXd> & coulomb_mat) {
        // initialize the distances to be sorted. the center is always first
        std::vector<double> distances_to_sort{};
        distances_to_sort.reserve(coulomb_mat.cols());

        auto row_norms = coulomb_mat.colwise().squaredNorm().eval();
        row_norms(0) = 1e200;
        for (auto idx_i{0}; idx_i < coulomb_mat.cols(); ++idx_i) {
          distances_to_sort.push_back(row_norms(idx_i));
        }

        std::vector<std::pair<size_t, distiter>> order_coulomb(
            distances_to_sort.size());
        size_t nn{0};
        for (distiter it{distances_to_sort.begin()};
             it != distances_to_sort.end(); ++it, ++nn) {
          order_coulomb[nn] = make_pair(nn, it);
        }

        // use stable sort
        std::stable_sort(order_coulomb.begin(), order_coulomb.end(),
                         ordering::descending);

        return order_coulomb;
      }
    };
    /* -------------------- rep-options-impl-end -------------------- */

  }  // namespace internal
  /* ---------------------------------------------------------------------- */
  /* -------------------- rep-preamble-start -------------------- */
  /**
   * Implementation of the Environmental Coulomb Matrix
   */
  template <class StructureManager>
  class RepresentationManagerSortedCoulomb : public RepresentationManagerBase {
   public:
    using Manager_t = StructureManager;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Parent = RepresentationManagerBase;
    // type of the hyperparameters
    using Hypers_t = typename Parent::Hypers_t;
    // numeric type for the representation features
    using Precision_t = typename Parent::Precision_t;
    // type of the data structure for the representation feaures
    using Property_t = Property<Precision_t, 1, 1,Manager_t, Eigen::Dynamic, 1>;
    template <size_t Order>
    // short hand type to help the iteration over the structure manager
    using ClusterRef_t = typename Manager_t::template ClusterRef<Order>;
    // type of the datastructure used to register the list of valid
    // hyperparameters
    using ReferenceHypers_t = Parent::ReferenceHypers_t;
    /* -------------------- rep-preamble-end -------------------- */

    /* -------------------- rep-construc-start -------------------- */
    //! Constructor
    RepresentationManagerSortedCoulomb(ManagerPtr_t sm, const Hypers_t & hyper)
        : structure_manager{std::move(sm)}, central_decay{},
          interaction_cutoff{}, interaction_decay{}, coulomb_matrices{*sm} {
      this->check_hyperparameters(this->reference_hypers, hyper);
      // Extract the options and hyperparameters
      this->set_hyperparameters(hyper);
      // additional checks specific to the coulomb matrix representation
      this->check_size_compatibility();
    }

    //! Copy constructor
    RepresentationManagerSortedCoulomb(
        const RepresentationManagerSortedCoulomb & other) = delete;

    //! Move constructor
    RepresentationManagerSortedCoulomb(
        RepresentationManagerSortedCoulomb && other) = default;

    //! Destructor
    virtual ~RepresentationManagerSortedCoulomb() = default;

    //! Copy assignment operator
    RepresentationManagerSortedCoulomb &
    operator=(const RepresentationManagerSortedCoulomb & other) = delete;

    //! Move assignment operator
    RepresentationManagerSortedCoulomb &
    operator=(RepresentationManagerSortedCoulomb && other) = default;
    /* -------------------- rep-construc-end -------------------- */

    /* -------------------- rep-interface-start -------------------- */
    //! compute representation
    void compute();

    //! set hypers
    void set_hyperparameters(const Hypers_t &);

    //! getter for the representation
    Eigen::Map<const Eigen::MatrixXd> get_representation_full() {
      auto nb_centers{this->structure_manager->size()};
      auto nb_features{this->get_n_feature()};
      auto & raw_data{this->coulomb_matrices.get_raw_data()};
      Eigen::Map<const Eigen::MatrixXd> representation(raw_data.data(),
                                                       nb_features, nb_centers);
      return representation;
    }

    //! get the raw data of the representation
    std::vector<Precision_t> & get_representation_raw_data() {
      return this->coulomb_matrices.get_raw_data();
    }

    Data_t & get_representation_sparse_raw_data() { return this->dummy; }

    //! get the size of a feature vector
    size_t get_feature_size() { return this->coulomb_matrices.get_nb_comp(); }

    //! get the number of centers for the representation
    size_t get_center_size() { return this->coulomb_matrices.get_nb_item(); }
    /* -------------------- rep-interface-end -------------------- */

    //! Implementation of compute representation
    template <internal::CMSortAlgorithm AlgorithmType>
    void compute_helper();

    //! check if size of representation manager is enough for current structure
    //! manager
    void check_size_compatibility() {
      for (auto center : this->structure_manager) {
        auto n_neighbours{center.size()};
        if (n_neighbours > this->size) {
          std::cout << "size is too small for this "
                       "structure and has been reset to: "
                    << n_neighbours << std::endl;
          this->size = n_neighbours;
        }
      }
    }

    //! returns the distance matrix for a central atom
    void get_distance_matrix(ClusterRef_t<1> & center,
                             Eigen::Ref<Eigen::MatrixXd> distance_mat,
                             Eigen::Ref<Eigen::MatrixXd> type_factor_mat);

    /**
     * Sort the coulomb matrix using the distance to the central atom
     * as reference order and linearize it.
     *
     * @params square_coulomb CM to sort
     * @params linear_coulomb sorted and linearized CM
     * @params order_coulomb  map to the sorted order order
     */
    template <typename DerivedA, typename DerivedB>
    void sort_and_linearize_coulomb_matrix(
        const Eigen::DenseBase<DerivedA> & square_coulomb,
        Eigen::DenseBase<DerivedB> & linear_coulomb,
        const std::vector<std::pair<size_t, internal::distiter>> &
            order_coulomb) {
      auto Nneigh{square_coulomb.cols()};
      size_t lin_id{0};

      for (int idx_i{0}; idx_i < Nneigh; ++idx_i) {
        size_t idx_is{order_coulomb[idx_i].first};
        for (int idx_j{0}; idx_j < idx_i + 1; ++idx_j) {
          size_t idx_js{order_coulomb[idx_j].first};
          linear_coulomb(lin_id) = square_coulomb(idx_is, idx_js);
          lin_id += 1;
        }
      }
    }

    //! scale cutoff factor depending on distance and decay
    inline double get_cutoff_factor(const double & distance,
                                    const double & cutoff,
                                    const double & decay) {
      if (distance <= cutoff - decay) {
        return 1.;
      } else if (distance > cutoff) {
        return 0.;
      } else if (distance > cutoff - decay) {
        return 0.5 *
               (1 + std::cos(EIGEN_PI * (distance - cutoff + decay) / decay));
      } else {
        throw std::runtime_error("Something went wrong...");
      }
    }

    //! get the size of a feature vector from the hyper parameters
    inline size_t get_n_feature() { return this->size * (this->size + 1) / 2; }

    /* -------------------- rep-variables-start -------------------- */
    // Reference to the structure manager
    ManagerPtr_t structure_manager;
    // list of hyperparameters specific to the coulomb matrix
    // spherical cutoff for the atomic environment
    double central_cutoff{};

    double central_decay{};
    double interaction_cutoff{};
    double interaction_decay{};
    // at least equal to the largest number of neighours
    size_t size{};

    Data_t dummy{};

    Property_t coulomb_matrices;

    //! reference the requiered hypers
    ReferenceHypers_t reference_hypers{
        {"central_decay", {}},
        {"interaction_cutoff", {}},
        {"interaction_decay", {}},
        {"size", {}},
        {"sorting_algorithm", {"distance", "row_norm"}},
    };
    /*-------------------- rep-variables-end -------------------- */
  };

  /* ---------------------------------------------------------------------- */
  template <class Mngr>
  void RepresentationManagerSortedCoulomb<Mngr>::set_hyperparameters(
      const RepresentationManagerSortedCoulomb<Mngr>::Hypers_t & hyper) {
    this->hypers = hyper;
    this->central_cutoff = this->structure_manager->get_cutoff();
    this->hypers["central_cutoff"] = this->central_cutoff;

    this->options.emplace("sorting_algorithm",
                          hyper["sorting_algorithm"].get<std::string>());

    this->size = hyper["size"];

    if ((hyper["interaction_cutoff"] < 0) or
        hyper["interaction_cutoff"] > 2 * this->central_cutoff) {
      this->interaction_cutoff = 2 * this->central_cutoff;
      this->hypers["interaction_cutoff"] = this->interaction_cutoff;
    } else {
      this->interaction_cutoff = hyper["interaction_cutoff"];
    }

    if (hyper["central_decay"] < 0) {
      this->central_decay = 0.;
      this->hypers["central_decay"] = this->central_decay;
    } else if (hyper["central_decay"] > this->central_cutoff) {
      this->central_decay = this->central_cutoff;
      this->hypers["central_decay"] = this->central_cutoff;
    } else {
      this->central_decay = hyper["central_decay"];
    }

    if (hyper["interaction_decay"] < 0) {
      this->interaction_decay = 0.;
      this->hypers["interaction_decay"] = this->interaction_decay;
    } else if (hyper["interaction_decay"] > this->interaction_cutoff) {
      this->interaction_decay = this->interaction_cutoff;
      this->hypers["interaction_decay"] = this->interaction_cutoff;
    } else {
      this->interaction_decay = hyper["interaction_decay"];
    }
  }

  /* ---------------------------------------------------------------------- */
  /* -------------------- rep-options-compute-start-------------------- */
  template <class Mngr>
  void RepresentationManagerSortedCoulomb<Mngr>::compute() {
    auto option{this->options["sorting_algorithm"]};

    if (option == "distance") {
      compute_helper<internal::CMSortAlgorithm::Distance>();
    } else if (option == "row_norm") {
      compute_helper<internal::CMSortAlgorithm::RowNorm>();
    } else {
      auto error_message{std::string("Option '") + option +
                         std::string("' is not implemented.")};
      throw std::invalid_argument(error_message.c_str());
    }
  }
  /* -------------------- rep-options-compute-end -------------------- */
  /* ---------------------------------------------------------------------- */
  /* -------------------- rep-options-compute-impl-start -------------------- */
  template <class Mngr>
  template <internal::CMSortAlgorithm AlgorithmType>
  void RepresentationManagerSortedCoulomb<Mngr>::compute_helper() {
    // initialise the sorted coulomb_matrices in linear storage
    this->coulomb_matrices.resize_to_zero();
    this->coulomb_matrices.set_nb_row(this->get_n_feature());

    // initialize the sorted linear coulomb matrix
    Eigen::MatrixXd lin_sorted_coulomb_mat(this->size * (this->size + 1) / 2,
                                           1);
    // Eigen::MatrixXd coulomb_mat(this->size,this->size);

    // loop over the centers
    for (auto center : *this->structure_manager) {
      // re-use the temporary coulomb mat in linear storage
      // need to be zeroed because old data might not be overwritten
      lin_sorted_coulomb_mat =
          Eigen::MatrixXd::Zero(this->size * (this->size + 1) / 2, 1);

      // n_neighbour counts the central atom and the neighbours
      size_t n_neighbour{center.size() + 1};

      // the local distance matrix. Ones to avoid overflow in the div.
      Eigen::MatrixXd distance_mat =
          Eigen::MatrixXd::Ones(n_neighbour, n_neighbour);
      // the matrix with the prefactor on the .
      Eigen::MatrixXd type_factor_mat =
          Eigen::MatrixXd::Zero(n_neighbour, n_neighbour);

      // the local coulomb matrix
      Eigen::MatrixXd coulomb_mat =
          Eigen::MatrixXd::Ones(n_neighbour, n_neighbour);

      this->get_distance_matrix(center, distance_mat, type_factor_mat);

      // Compute Coulomb Mat element wise.
      coulomb_mat = type_factor_mat.array() / distance_mat.array();

      using Sorter = internal::SortCoulomMatrix<AlgorithmType>;
      auto sort_order{
          Sorter::get_coulomb_matrix_sorting_order(distance_mat, coulomb_mat)};
      // inject the coulomb matrix into the sorted linear storage
      this->sort_and_linearize_coulomb_matrix(
          coulomb_mat, lin_sorted_coulomb_mat, sort_order);

      this->coulomb_matrices.push_back(lin_sorted_coulomb_mat);
    }
  }
  /* -------------------- rep-options-compute-impl-end -------------------- */

  /* ---------------------------------------------------------------------- */
  template <class Mngr>
  void RepresentationManagerSortedCoulomb<Mngr>::get_distance_matrix(
      RepresentationManagerSortedCoulomb<Mngr>::ClusterRef_t<1> & center,
      Eigen::Ref<Eigen::MatrixXd> distance_mat,
      Eigen::Ref<Eigen::MatrixXd> type_factor_mat) {
    // the coulomb mat first row and col corresponds
    // to central atom to neighbours
    auto && Zk{center.get_atom_type()};
    auto && central_cutoff{this->structure_manager->get_cutoff()};

    type_factor_mat(0, 0) = 0.5 * std::pow(Zk, 2.4);
    for (auto neigh_i : center) {
      size_t idx_i{neigh_i.get_index() + 1};
      auto && Zi{neigh_i.get_atom_type()};
      double & dik{this->structure_manager->get_distance(neigh_i)};
      double fac_ik{
          get_cutoff_factor(dik, central_cutoff, this->central_decay)};

      type_factor_mat(idx_i, 0) = Zk * Zi * fac_ik * fac_ik;
      type_factor_mat(idx_i, idx_i) = 0.5 * std::pow(Zi, 2.4) * fac_ik * fac_ik;
      distance_mat(idx_i, 0) = dik;
      type_factor_mat(0, idx_i) = type_factor_mat(idx_i, 0);
      distance_mat(0, idx_i) = distance_mat(idx_i, 0);
    }

    // compute the neighbour to neighbour part of the coulomb matrix
    for (auto neigh_i : center) {
      size_t idx_i{neigh_i.get_index() + 1};
      auto && Zi{neigh_i.get_atom_type()};
      double & dik{this->structure_manager->get_distance(neigh_i)};
      double fac_ik{
          get_cutoff_factor(dik, central_cutoff, this->central_decay)};

      for (auto neigh_j : center) {
        size_t idx_j{neigh_j.get_index() + 1};
        // work only on the lower diagonal
        if (idx_i >= idx_j)
          continue;
        auto && Zj{neigh_j.get_atom_type()};
        double dij{(neigh_i.get_position() - neigh_j.get_position()).norm()};
        double fac_ij{get_cutoff_factor(dij, this->interaction_cutoff,
                                        this->interaction_decay)};

        double djk{(center.get_position() - neigh_j.get_position()).norm()};
        double fac_jk{get_cutoff_factor(djk, this->interaction_cutoff,
                                        this->interaction_decay)};

        type_factor_mat(idx_j, idx_i) = Zj * Zi * fac_ij * fac_ik * fac_jk;
        distance_mat(idx_j, idx_i) = dij;
        type_factor_mat(idx_i, idx_j) = type_factor_mat(idx_j, idx_i);
        distance_mat(idx_i, idx_j) = distance_mat(idx_j, idx_i);
      }
    }
  }
}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SORTED_COULOMB_HH_
