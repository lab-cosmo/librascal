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

  template <SymmetryFunType SymFunType, class StructureManager>
  void BehlerFeatureBase::compute_helper(StructureManager & manager) const {
    switch (this->cut_fun_type) {
    case CutoffFunctionType::Cosine: {
      auto & feature{static_cast<
          const BehlerFeature<SymFunType, CutoffFunctionType::Cosine> &>(
          *this)};
      feature.compute(manager);
      break;
    }
    case CutoffFunctionType::RadialScaling: {
      auto & feature{static_cast<
          const BehlerFeature<SymFunType, CutoffFunctionType::RadialScaling> &>(
          *this)};
      feature.compute(manager);
      break;
    }
    default:
      throw std::runtime_error("Unknown cutoff function type");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  void BehlerFeatureBase::compute(StructureManager & manager) const {
    switch (this->sym_fun_type) {
    case SymmetryFunType::One: {
      this->compute_helper<SymmetryFunType::One>(manager);
      break;
    }
    case SymmetryFunType::Gaussian: {
      this->compute_helper<SymmetryFunType::Gaussian>(manager);
      break;
    }
    case SymmetryFunType::Cosine: {
      this->compute_helper<SymmetryFunType::Cosine>(manager);
      break;
    }
    case SymmetryFunType::Angular1: {
      this->compute_helper<SymmetryFunType::Angular1>(manager);
      break;
    }
    case SymmetryFunType::Angular2: {
      this->compute_helper<SymmetryFunType::Angular2>(manager);
      break;
    }
    default:
      throw std::runtime_error("Unknown symmetry function type");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, internal::CutoffFunctionType CutFunType>
  template <class StructureManager>
  void BehlerFeature<SymFunType, CutFunType>::compute(
      StructureManager & /*manager*/) const {
    // utils::for_each_at_order<SymmetryFun<SymFunType>::NbParams>::loop(
    //     eval_cluster, manager);
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, internal::CutoffFunctionType CutFunType>
  void BehlerFeature<SymFunType, CutFunType>::init(const UnitStyle & units) {
    std::map<double, size_t> nb_param_per_cutoff{};

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
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_IMPL_HH_
