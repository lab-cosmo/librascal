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
  using internal::CutoffFunction;
  using internal::CutoffFunctionType;

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType... SymFunTypes>
  class BehlerFeatureBase::SymFunctionsVTable {};

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType Head, SymmetryFunType... Tail>
  class BehlerFeatureBase::SymFunctionsVTable<Head, Tail...> {
    template <class StructureManager, class... PropertyPtr>
    static void compute(BehlerFeatureBase & behler_feature,
                        StructureManager & manager, PropertyPtr... outputs) {
      if (behler_feature.sym_fun_type == Head) {
        behler_feature.compute_helper<Head>(manager, outputs...);
      } else {
        SymFunctionsVTable<Tail...>::compute(behler_feature, manager,
                                             outputs...);
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType Head>
  class BehlerFeatureBase::SymFunctionsVTable<Head> {
    template <class StructureManager, class... PropertyPtr>
    static void compute(BehlerFeatureBase & behler_feature,
                        StructureManager & manager, PropertyPtr... outputs) {
      if (behler_feature.sym_fun_type == Head) {
        behler_feature.compute_helper<Head>(manager, outputs...);
      } else {
        std::stringstream err{};
        err << "Symmetry function type " << behler_feature.sym_fun_type
            << " is not known.";
        throw std::runtime_error(err.str());
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, CutoffFunctionType... CutoffFunTypes>
  class BehlerFeatureBase::CutoffFunctionsVTable {};

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, CutoffFunctionType Head,
            CutoffFunctionType... Tail>
  class BehlerFeatureBase::CutoffFunctionsVTable<SymFunType, Head, Tail...> {
    template <class StructureManager, class... PropertyPtr>
    static void compute(BehlerFeatureBase & behler_feature,
                        StructureManager & manager, PropertyPtr... outputs) {
      if (behler_feature.cut_fun_type == Head) {
        auto & feature{dynamic_cast<
            const BehlerFeature<SymFunType, CutoffFunctionType::Cosine> &>(
            *behler_feature)};
        feature.compute(manager, outputs...);
      } else {
        CutoffFunctionsVTable<SymFunType, Tail...>::compute(
            behler_feature, manager, outputs...);
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, CutoffFunctionType Head>
  class BehlerFeatureBase::CutoffFunctionsVTable<SymFunType, Head> {
    if (behler_feature.cut_fun_type == Head) {
      auto & feature{dynamic_cast<
          const BehlerFeature<SymFunType, CutoffFunctionType::Cosine> &>(
          *behler_feature)};
      feature.compute(manager, outputs...);
    } else {
      std::stringstream err{};
      err << "Cutoff function type " << behler_feature.cut_fun_type
          << " is not known.";
      throw std::runtime_error(err.str());
    }
  };

  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  void BehlerFeatureBase::compute(StructureManager & manager,
                                  std::shared_ptr<PropertyBase> prop) const {
    SymFunctionsVTable<SymmetryFunType::One, SymmetryFunType::Gaussian,
                       SymmetryFunType::Cosine, SymmetryFunType::Angular1,
                       SymmetryFunType::Angular2>::compute(*this, manager,
                                                           prop);
  }

  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  void
  BehlerFeatureBase::compute(StructureManager & manager,
                             std::shared_ptr<PropertyBase> prop,
                             std::shared_ptr<PropertyBase> prop_der) const {
    SymFunctionsVTable<SymmetryFunType::One, SymmetryFunType::Gaussian,
                       SymmetryFunType::Cosine, SymmetryFunType::Angular1,
                       SymmetryFunType::Angular2>::compute(*this, manager, prop,
                                                           prop_der);
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, class StructureManager>
  void BehlerFeatureBase::compute_helper(StructureManager & manager,
                                         std::shared_ptr<PropertyBase> output) {
    CutoffFunctionsVTable<SymFunType, SymmetryFunType::One,
                          SymmetryFunType::Gaussian, SymmetryFunType::Cosine,
                          SymmetryFunType::Angular1,
                          SymmetryFunType::Angular2>::compute(*this, manager,
                                                              output);
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, class StructureManager>
  void BehlerFeatureBase::compute_helper(
      StructureManager & manager, std::shared_ptr<PropertyBase> output,
      std::shared_ptr<PropertyBase> output_derivatives) {
    CutoffFunctionsVTable<
        SymFunType, SymmetryFunType::One, SymmetryFunType::Gaussian,
        SymmetryFunType::Cosine, SymmetryFunType::Angular1,
        SymmetryFunType::Angular2>::compute(*this, manager, output,
                                            output_derivatives);
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, internal::CutoffFunctionType CutFunType>
  template <class StructureManager>
  void BehlerFeature<SymFunType, CutFunType>::compute(
      StructureManager & manager, std::shared_ptr<PropertyBase> output) const {
    const std::string && cut_off_vals_id{this->cut_off_fun.get_identifier()};
    auto & cutoffs{
        dynamic_cast<Property<> &>(*manager.get_property(cut_off_vals_id))};
    if (not manager.get_property
    using Output_t =
        Property<double, AtomOrder, PropertyLayer, StructureManager>;
    Output_t & fun_vals{dynamic_cast<Output_t &>(*output)};
    for (auto && atom : manager) {
      for (auto && cluster :
           atom.template get_clusters_of_order<SymmetryFunction::Order>()) {
        this->eval_cluster(cluster, fun_vals[atom])
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, internal::CutoffFunctionType CutFunType>
  template <class StructureManager>
  void BehlerFeature<SymFunType, CutFunType>::compute(
      StructureManager & manager, std::shared_ptr<PropertyBase> output,
      std::shared_ptr<PropertyBase> output) const {
    using Output_t = Property<double, SymmetryFunction::Order, PropertyLayer,
                              StructureManager>;
    Output_t & fun_vals{dynamic_cast<Output_t &>(*output)};
    if (not_equal) {
      for (auto && atom : manager) {
        for (auto && cluster :
             atom.template get_clusters_of_order<SymmetryFunction::Order>()) {
          fun_vals[atom] = this->eval_cluster(cluster)
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, internal::CutoffFunctionType CutFunType>
  void BehlerFeature<SymFunType, CutFunType>::init(const UnitStyle & units) {
    if (this->is_initialised) {
      throw std::runtime_error("double initialisation");
    }
    // counting the number of parameters to store per cutoff
    for (const auto & param : this->raw_params) {
      auto && r_cut{param.at("r_cut").template get<double>()};
      nb_param_per_cutoff[r_cut]++;
    }

    // allocate storage
    for (const auto & key_val : nb_param_per_cutoff) {
      auto && r_cut{key_val.first};
      this->params[r_cut].resize(SymmetryFun<SymFunType>::NbParams,
                                 nb_param_per_cutoff[r_cut]);
      this->params[r_cut].setZero();
    }

    // store params in storage
    std::map<double, size_t> nb_param_counter{};
    for (auto && param : this->raw_params) {
      auto && r_cut{param.at("r_cut").template get<double>()};

      this->params[r_cut].col(nb_param_counter[r_cut]++) =
          SymmetryFun<SymFunType>::read(param, units);
    }

    // determine which species repetition scenario we are in
    for (auto && species : param.at("species").template get<std::string>()) {
      this->species_combo
    }
    auto species{};

    this->is_initialised = true;
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_IMPL_HH_
