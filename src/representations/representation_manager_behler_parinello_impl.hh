/**
 * file   representation_manager_behler_parinello_impl.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   13 Dec 2018
 *
 * @brief  implementation for Behler-Parinello representation manager
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with librascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "utils/for_each_at_order.hh"
namespace rascal {

  template <class StructureManager>
  BehlerParinello<StructureManager>::BehlerParinello(
      StructureManager & structure, const json & hypers)
      : structure{structure}, species{structure} {}

  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  void BehlerParinello<StructureManager>::update() {
    this->species.update();
    for (auto && tup : this->species.filters_by_nb_elements<1>()) {
      auto && species{tup.first};
      auto && filter{tup.second};

      // make sure all storage properties exist
      using GProperty_t = TypedProperty<double, AtomOrder>;
      std::shared_ptr<GProperty_t> G_values{}, dG_values{};
      if (not filter.has_property(this->symmetry_function_key)) {
        constexpr size_t AtomOrder{1};

        G_values = std::make_shared<GProperty_t>(
            filter, this->nb_sym_per_species[std::get<0>(species)]);
        filter.attach_property(this->symmetry_function_key, G_values);

        dG_values = std::make_shared<GProperty_t>(
            filter, this->nb_sym_per_species[std::get<0>(species)], Dim);
        filter.attach_property(this->symmetry_derivative_key, G_values);
      } else {
        G_values = filter.get_property(this->symmetry_function_key);
        dG_values = filter.get_property(this->symmetry_derivative_key);
      }

      // resize storage properties to updated species managers
      G_values.resize();
      dG_values.resize();
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  void BehlerParinello<StructureManager>::compute() {
    using utils::for_each_at_order;

    //! precompute precomputable values
    for (const auto && species_key_val: this->symmetry_functions) {
      auto && sym_fun{species_key_val.second};
      sym_fun.prepare(this->structure);
    }

    constexpr auto PairOrder{2};
    constexpr auto TripletOrder{3};
    for (const auto && species_key_val : this->symmetry_functions) {
      auto && species_combo{species_key_val.first};
      auto && sym_fun{species_key_val.second};

      const auto & order{species_combo.get_order()};

      switch (order) {
      case PairOrder: {
        this->evaluate_pair_symmetry_function_group(sym_fun);
        break;
      }
      case TripletOrder: {
        this->evaluate_triplet_symmetry_function_group(sym_fun);
        break;
      }
      default:
        break;
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  template <size_t Order>
  void BehlerParinello<StructureManager>::evaluate_pair_symmetry_function_group(
      std::vector<InputNodeContributionBase> & symmetry_funs) {
    for (auto && symmetry_fun : symmetry_funs) {
      this->evaluate_sym_fun<Permutation<0, 1>>((symmetry_fun));
      this->evaluate_sym_fun<Permutation<1, 0>>((symmetry_fun));
    }
  }
  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  template <size_t Order>
  void
  BehlerParinello<StructureManager>::evaluate_triplet_symmetry_function_group(
      std::vector<InputNodeContributionBase> & symmetry_funs) {
    for (auto && symmetry_fun : symmetry_funs) {
      this->evaluate_sym_fun<Permutation<0, 1, 2>>((symmetry_fun));
      this->evaluate_sym_fun<Permutation<1, 2, 0>>((symmetry_fun));
      this->evaluate_sym_fun<Permutation<2, 0, 1>>((symmetry_fun));

      if (not this->legacy_behaviour) {
        this->evaluate_sym_fun<Permutation<2, 1, 0>>((symmetry_fun));
        this->evaluate_sym_fun<Permutation<1, 0, 2>>((symmetry_fun));
        this->evaluate_sym_fun<Permutation<0, 2, 1>>((symmetry_fun));
      }
    }
  }

}  // namespace rascal
