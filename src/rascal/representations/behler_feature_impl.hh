/**
 * file   behler_feature_impl.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief  implementation for input node contribution tools
 *
 * Copyright © 2019 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_IMPL_HH_
#define SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_IMPL_HH_

#include "rascal/utils/for_each_at_order.hh"

namespace rascal {
  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType... SymFunTypes>
  template <SymmetryFunctionType... SymFunTypes_>
  class BehlerFeatureBase<SymFunTypes...>::SymFunctionsVTable {};

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType... SymFunTypes>
  template <SymmetryFunctionType Head, SymmetryFunctionType... Tail>
  class BehlerFeatureBase<SymFunTypes...>::SymFunctionsVTable<Head, Tail...> {
    template <class StructureManager, class... PropertyPtr>
    static void compute(BehlerFeatureBase & behler_feature,
                        StructureManager & manager, PropertyPtr... outputs) {
      if (behler_feature.sym_fun_type == Head) {
        auto & feature{
            dynamic_cast<const BehlerFeature<Head, SymFunTypes...> &>(
                behler_feature)};
        feature.compute(manager, outputs...);
      } else {
        SymFunctionsVTable<Tail...>::compute(behler_feature, manager,
                                             outputs...);
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType... SymFunTypes>
  template <SymmetryFunctionType Head>
  class BehlerFeatureBase<SymFunTypes...>::SymFunctionsVTable<Head> {
    template <class StructureManager, class... PropertyPtr>
    static void compute(BehlerFeatureBase & behler_feature,
                        StructureManager & manager, PropertyPtr... outputs) {
      if (behler_feature.sym_fun_type == Head) {
        auto & feature{
            dynamic_cast<const BehlerFeature<Head, SymFunTypes...> &>(
                behler_feature)};
        feature.compute(manager, outputs...);
      } else {
        std::stringstream err{};
        err << "Symmetry function type " << behler_feature.sym_fun_type
            << " is not known.";
        throw std::runtime_error(err.str());
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType... SymFunTypes>
  template <class StructureManager>
  void BehlerFeatureBase<SymFunTypes...>::compute(
      StructureManager & manager, std::shared_ptr<PropertyBase> prop) const {
    SymFunctionsVTable<
        SymmetryFunctionType::One, SymmetryFunctionType::Gaussian,
        SymmetryFunctionType::Cosine, SymmetryFunctionType::Angular1,
        SymmetryFunctionType::Angular2>::compute(*this, manager, prop);
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType... SymFunTypes>
  template <class StructureManager>
  void BehlerFeatureBase<SymFunTypes...>::compute(
      StructureManager & manager, std::shared_ptr<PropertyBase> prop,
      std::shared_ptr<PropertyBase> prop_der) const {
    SymFunctionsVTable<
        SymmetryFunctionType::One, SymmetryFunctionType::Gaussian,
        SymmetryFunctionType::Cosine, SymmetryFunctionType::Angular1,
        SymmetryFunctionType::Angular2>::compute(*this, manager, prop,
                                                 prop_der);
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  template <class StructureManager, RepeatedSpecies RepSpecies,
            typename Permutation>
  void BehlerFeature<MySymFunType, SymFunTypes...>::compute_helper(
      StructureManager & manager, std::shared_ptr<PropertyBase> output) const {
    const std::string && cut_off_vals_id{this->cut_off_fun.get_identifier()};
    auto & cutoffs{this->cut_fun.get_value(manager)};
    // eval
    using Output_t = Property<double, AtomOrder, StructureManager>;
    Output_t & fun_vals{dynamic_cast<Output_t &>(*output)};
    auto & distances{manager->get_distance()};

    switch (SymmetryFunction_t::Order) {
    case PairOrder: {
      for (auto && atom : manager) {
        for (auto && pair : atom.pairs()) {
          auto && leading_cluster_id{Permutation::leading(manager, pair)};
          auto && G_incr{this->sym_fun.f_sym(distances[pair]) * cutoffs[pair]};
          if (Permutation::no) {
            fun_vals[leading_cluster_id] += G_incr;
          } else {
            fun_vals[pair] += G_incr;
          }
          switch (RepSpecies) {
          case RepeatedSpecies::Not: {
            break;
          }
          case RepeatedSpecies::All: {
            auto && second_cluster_id{Permutation::second(manager, pair)};
            fun_vals[second_cluster_id] += G_incr;
            break;
          }
          default:
            throw std::runtime_error("Unknown species repetition pattern");
            break;
          }
        }
      }
      break;
    }
    case TripletOrder: {
      break;
    }
    default:
      std::stringstream err{};
      err << "unknown symmetry function order " << SymmetryFunction_t::Order;
      throw std::runtime_error(err.str());
      break;
    }
  }  // namespace rascal
  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  template <RepeatedSpecies RepSpecies, class StructureManager>
  void BehlerFeature<MySymFunType, SymFunTypes...>::compute(
      StructureManager & manager, std::shared_ptr<PropertyBase> output_value,
      std::shared_ptr<PropertyBase> output_derivative) const {
    using Output_t =
        Property<double, SymmetryFunction_t::Order, StructureManager>;
    auto & cutoff_values{this->cut_fun.get_value(manager)};
    auto & cutoff_derivatives{this->cut_fun.get_derivative(manager)};

    Output_t & fun_values{dynamic_cast<Output_t &>(*output_value)};
    Output_t & fun_derivatives{dynamic_cast<Output_t &>(*output_derivative)};
    for (auto && atom : manager) {
      for (auto && cluster :
           atom.template get_clusters_of_order<SymmetryFunction_t::Order>()) {
        std::tie(fun_values[atom], fun_derivatives[cluster])  = this->eval_cluster(cluster);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  void
  BehlerFeature<MySymFunType, SymFunTypes...>::init(const UnitStyle & /*units*/) {
    // if (this->is_initialised) {
    //   throw std::runtime_error("double initialisation");
    // }
    // // counting the number of parameters to store per cutoff
    // for (const auto & param : this->raw_params) {
    //   auto && r_cut{param.at("r_cut").template get<double>()};
    //   nb_param_per_cutoff[r_cut]++;
    // }

    // // allocate storage
    // for (const auto & key_val : nb_param_per_cutoff) {
    //   auto && r_cut{key_val.first};
    //   this->params[r_cut].resize(SymmetryFun<SymFunType>::NbParams,
    //                              nb_param_per_cutoff[r_cut]);
    //   this->params[r_cut].setZero();
    // }

    // // store params in storage
    // std::map<double, size_t> nb_param_counter{};
    // for (auto && param : this->raw_params) {
    //   auto && r_cut{param.at("r_cut").template get<double>()};

    //   this->params[r_cut].col(nb_param_counter[r_cut]++) =
    //       SymmetryFun<SymFunType>::read(param, units);
    // }

    // // determine which species repetition scenario we are in
    // for (auto && species : param.at("species").template get<std::string>()) {
    //   this->species_combo
    // }
    // auto species{};

    // this->is_initialised = true;
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_IMPL_HH_
