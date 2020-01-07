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

  template <class Calculator>
  class PseudoPointsBlockSparse {
   public:
    using Key_t = typename CalculatorBase::Key_t;
    using Data_t = std::map<int, std::map<Key_t, std::vector<double>>>;
    using Indices_t = std::map<int, std::map<Key_t, std::vector<size_t>>>;
    using Counters_t = std::map<int, size_t>;

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

    PseudoPointsBlockSparse() {
      // there is less than 130 elemtents
      for (int sp{1}; sp < 130; ++sp) {
        counters[sp] = 0;
      }
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
