/**
 * @file   rascal/models/pseudo_points.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   18 Dec 2019
 *
 * @brief Implementation of pseudo points for sparse kernels
 *
 * Copyright 2019 Felix Musil COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_RASCAL_MODELS_PSEUDO_POINTS_HH_
#define SRC_RASCAL_MODELS_PSEUDO_POINTS_HH_

#include "rascal/math/utils.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/structure_managers/property_block_sparse.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/utils/json_io.hh"

namespace rascal {

  namespace internal {}  // namespace internal

  /**
   * Set of pseudo points associated with a representation stored in BlockSparse
   * format, e.g. SphericalInvariants. The number of pseudo points is often
   * refered as $M$ and note that they might be the representation of actual
   * atomic environments or completely artificial.
   */
  template <class Calculator>
  class PseudoPointsBlockSparse {
   public:
    using Key_t = typename CalculatorBase::Key_t;
    using Data_t = std::map<int, std::map<Key_t, std::vector<double>>>;
    using Indices_t = std::map<int, std::map<Key_t, std::vector<size_t>>>;
    using Counters_t = std::map<int, size_t>;
    using ColVector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using ColVectorDer_t = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::ColMajor>;

    template <class StructureManager>
    using Property_t =
        typename Calculator::template Property_t<StructureManager>;
    template <class StructureManager>
    using PropertyGradient_t =
        typename Calculator::template PropertyGradient_t<StructureManager>;

    Data_t values{};
    Indices_t indicies{};
    Counters_t counters{};
    size_t inner_size{0};
    std::set<int> center_species{};
    std::set<Key_t> keys{};

    constexpr static int n_spatial_dimensions{3};

    PseudoPointsBlockSparse() {
      // there is less than 130 elemtents
      for (int sp{1}; sp < 130; ++sp) {
        counters[sp] = 0;
      }
    }

    size_t size() const {
      size_t n_points{0};
      for (const auto& count: this->counters) {
        n_points += count.second;
      }
      return n_points;
    }

    const std::set<int> & species() const {
      return center_species;
    }

    std::vector<int> species_by_points() const {
      std::vector<int> species{};
      for (const auto & sp : this->center_species) {
        for (size_t ii{0}; ii < counters.at(sp); ++ii) {
          species.push_back(sp);
        }
      }
      return species;
    }

    math::Matrix_t self_dot(const int& sp) const {
      math::Matrix_t KMM_by_sp(this->counters.at(sp), this->counters.at(sp));
      KMM_by_sp.setZero();
      const auto & values_by_sp = values.at(sp);
      const auto & indicies_by_sp = indicies.at(sp);
      for (const Key_t& key : keys) {
        const auto & indicies_by_sp_key = indicies_by_sp.at(key);
        auto mat = Eigen::Map<const math::Matrix_t>(values_by_sp.at(key).data(),          static_cast<Eigen::Index>(indicies_by_sp_key.size()),
              static_cast<Eigen::Index>(this->inner_size));
        // auto KMM_by_key = (mat * mat.transpose()).eval();
        math::Matrix_t KMM_by_key(indicies_by_sp_key.size(), indicies_by_sp_key.size());
        KMM_by_key = KMM_by_key.setZero().selfadjointView<Eigen::Upper>().rankUpdate(mat);
        for (int i_row{0}; i_row < KMM_by_key.rows(); i_row++) {
          for (int i_col{0}; i_col < KMM_by_key.cols(); i_col++) {
            KMM_by_sp(indicies_by_sp_key[i_row], indicies_by_sp_key[i_col]) +=
                                                  KMM_by_key(i_row, i_col);
          }
        }
      } // key
      return KMM_by_sp;
    }

    /**
     * Compute the dot product between the pseudo points associated with type
     * sp with the representation associated with a center.
     * @return column vector Mx1
     */
    template <class V>
    ColVector_t dot(const int& sp, internal::InternallySortedKeyMap<Key_t, V> & representation) const {
      const auto & values_by_sp = this->values.at(sp);
      const auto & indicies_by_sp = this->indicies.at(sp);
      ColVector_t KNM_row(this->size());
      KNM_row.setZero();
      if (this->center_species.count(sp) == 0) {
        // the type of the central atom is not in the pseudo points
        return KNM_row;
      }
      int offset{0};
      for (const int& csp : this->center_species) {
        if (csp == sp) {
          break;
        } else {
          offset += this->counters.at(csp);
        }
      }

      for (const Key_t& key : this->keys) {
        if (representation.count(key)) {
          auto rep_flat_by_key{representation.flattened(key)};
          const auto & indicies_by_sp_key = indicies_by_sp.at(key);
          auto mat = Eigen::Map<const math::Matrix_t>(values_by_sp.at(key).data(),          static_cast<Eigen::Index>(indicies_by_sp_key.size()),
                static_cast<Eigen::Index>(this->inner_size));
          auto KNM_row_key = (mat * rep_flat_by_key.transpose()).eval();
          for (int i_row{0}; i_row < KNM_row_key.size(); i_row++) {
            KNM_row(offset+indicies_by_sp_key[i_row]) += KNM_row_key(i_row);
          }
        }
      }
      return KNM_row;
    }

    /**
     * Compute the dot product between the pseudo points associated with type
     * sp with the gradient of the representation associated with a center.
     * @return column vector Mx3
     */
    template <class V>
    ColVectorDer_t dot_derivative(const int& sp, internal::InternallySortedKeyMap<Key_t, V> & representation_grad) const {
      const auto & values_by_sp = this->values.at(sp);
      const auto & indicies_by_sp = this->indicies.at(sp);
      ColVectorDer_t KNM_row(this->size(), n_spatial_dimensions);
      KNM_row.setZero();
      if (this->center_species.count(sp) == 0) {
        // the type of the central atom is not in the pseudo points
        return KNM_row;
      }
      int offset{0};
      for (const int& csp : this->center_species) {
        if (csp == sp) {
          break;
        } else {
          offset += this->counters.at(csp);
        }
      }

      for (const Key_t& key : this->keys) {
        if (representation_grad.count(key)) {
          auto rep_grad_flat_by_key{representation_grad.flattened(key)};
          Eigen::Map<Eigen::Matrix<double, n_spatial_dimensions,
            Eigen::Dynamic, Eigen::RowMajor>> rep_grad_by_key(rep_grad_flat_by_key.data(), n_spatial_dimensions, this->inner_size);
          const auto & indicies_by_sp_key = indicies_by_sp.at(key);
          auto mat = Eigen::Map<const math::Matrix_t>(values_by_sp.at(key).data(),          static_cast<Eigen::Index>(indicies_by_sp_key.size()),
                static_cast<Eigen::Index>(this->inner_size));
          auto KNM_row_key = (mat * rep_grad_flat_by_key.transpose()).eval();
          for (int i_row{0}; i_row < KNM_row_key.rows(); i_row++) {
            KNM_row.row(offset+indicies_by_sp_key[i_row]) += KNM_row_key.row(i_row);
          }
        }
      }
      return KNM_row;
    }

    template <class ManagerCollection>
    void
    push_back(const Calculator & calculator,
              const ManagerCollection & collection,
              const std::vector<std::vector<int>> & selected_center_indices) {
      for (size_t i_manager{0}; i_manager < collection.size(); ++i_manager) {
        this->push_back(calculator, collection[i_manager],
                        selected_center_indices[i_manager]);
      }
    }

    template <class StructureManager>
    void push_back(const Calculator & calculator,
                   std::shared_ptr<StructureManager> manager,
                   const std::vector<int> & selected_center_indices) {
      std::string prop_name{calculator.get_name()};
      auto property =
          manager->template get_property<Property_t<StructureManager>>(
              prop_name);
      for (const auto & center_index : selected_center_indices) {
        auto center_it = manager->get_iterator_at(center_index);
        auto center = *center_it;
        auto && map2row = property->operator[](center);
        this->push_back(map2row, center.get_atom_type());
      }
    }

    template <class V>
    void push_back(internal::InternallySortedKeyMap<Key_t, V> & pseudo_point,
                   const int & center_type) {
      auto & values_by_sp = this->values[center_type];
      auto & counters_by_sp = this->counters[center_type];
      auto & indicies_by_sp = this->indicies[center_type];
      this->center_species.insert(center_type);
      for (const auto & key : pseudo_point.get_keys()) {
        this->keys.insert(key);
        auto & values_by_sp_key = values_by_sp[key];
        // Eigen::Map<const V> vals{pseudo_point[key]};
        // auto pseudo_point_by_key = Eigen::Map<math::Vector_t>(&vals(0, 0),
        // vals.size());
        auto pseudo_point_by_key = pseudo_point.flattened(key);
        for (int ii{0}; ii < pseudo_point_by_key.size(); ++ii) {
          values_by_sp_key.push_back(pseudo_point_by_key[ii]);
        }
        indicies_by_sp[key].push_back(counters_by_sp);
        if (this->inner_size == 0) {
          this->inner_size = pseudo_point_by_key.size();
        } else if (static_cast<Eigen::Index>(this->inner_size) !=
                   pseudo_point_by_key.size()) {
          std::stringstream err_str{};
          err_str << "The representation changed size during the set-up of "
                     "PseudoPointsBlockSparse:"
                  << "'" << this->inner_size
                  << "!=" << pseudo_point_by_key.size() << "'.";
          throw std::logic_error(err_str.str());
        }
      }
      ++counters_by_sp;
    }

    math::Matrix_t get_features() {
      size_t n_pseudo_points{0};
      for (const auto & sp : this->center_species) {
        n_pseudo_points += counters[sp];
      }
      math::Matrix_t mat{n_pseudo_points, this->inner_size * this->keys.size()};
      mat.setZero();
      size_t i_row{0};
      for (const auto & sp : this->center_species) {
        auto & values_by_sp = this->values[sp];
        auto & indicies_by_sp = this->indicies[sp];
        size_t i_col{0};
        for (const auto & key : this->keys) {
          if (values_by_sp.count(key)) {
            Eigen::Map<math::Matrix_t> block{
                values_by_sp[key].data(),
                static_cast<Eigen::Index>(indicies_by_sp[key].size()),
                static_cast<Eigen::Index>(this->inner_size)};
            for (size_t ii{0}; ii < indicies_by_sp[key].size(); ii++) {
              mat.block(i_row + indicies_by_sp[key][ii], i_col, 1,
                        this->inner_size) = block.row(ii);
            }
          }
          i_col += this->inner_size;
        }
        i_row += this->counters[sp];
      }
      return mat;
    }

    // math::Matrix_t
  };

}  // namespace rascal

#endif  // SRC_RASCAL_MODELS_PSEUDO_POINTS_HH_
