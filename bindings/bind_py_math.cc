/**
 * @file   bind_py_math.cc
 *
 * @author Michele Ceriotti <michele.ceriotti@gmail.com>
 *
 * @date   22 August 2018
 *
 * @brief
 *
 * Copyright  2018  Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
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

#include "bind_py_math.hh"

namespace rascal {

  void bind_sph(py::module & mod) {
    using math::SphericalHarmonics;
    py::class_<SphericalHarmonics> sph(mod, "SphericalHarmonics");
    sph.def(py::init([](size_t max_angular, bool calculate_derivatives) {
      auto sph_ = std::make_unique<SphericalHarmonics>(calculate_derivatives);
      sph_->precompute(max_angular);
      return sph_;
    }));
    sph.def("compute",
            py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(
                &SphericalHarmonics::calc));
    sph.def("harmonics", &SphericalHarmonics::get_harmonics,
            "return vector([l_max+1]**2)");
    sph.def("harmonics_derivatives",
            &SphericalHarmonics::get_harmonics_derivatives,
            "return matrix(3, [l_max+1]**2)");
  }

  void bind_ri(py::module & mod) {
    using RadialContribution_t =
        internal::RadialContribution<internal::RadialBasisType::DVR>;

    py::class_<RadialContribution_t> ri_dvr(mod, "RadialContributionDVR");

    ri_dvr.def(py::init([](const py::dict & hyper) {
      // convert to json
      json hypers = hyper;
      return std::make_unique<RadialContribution_t>(hypers);
    }));

    ri_dvr.def(
        "compute_center_contribution",
        [](RadialContribution_t & ri, const double fac_a) {
          return ri.compute_center_contribution(fac_a);
        },
        "fac_a = 1 / (2*sigma^2) and returns a vector(n_max)");

    ri_dvr.def(
        "compute_neighbour_contribution",
        [](RadialContribution_t & ri, const double distance,
           const double fac_a) {
          return ri.compute_neighbour_contribution(distance, fac_a);
        },
        "fac_a = 1 / (2*sigma^2) and returns a matrix(n_max, l_max+1)");

    ri_dvr.def(
        "compute_neighbour_derivative",
        [](RadialContribution_t & ri) {
          return ri.compute_neighbour_derivative();
        },
        "returns a matrix(n_max, l_max+1)");
  }

  void bind_fc(py::module & mod) {
    using CutoffFunction_t =
        internal::CutoffFunction<internal::CutoffFunctionType::ShiftedCosine>;

    py::class_<CutoffFunction_t> fc_c(mod, "CutoffFunctionShiftedCosine");
    fc_c.def(py::init([](const py::dict & hyper) {
               // convert to json
               json hypers = hyper;
               return std::make_unique<CutoffFunction_t>(hypers);
             }),
             "hyper should be a dict(cutoff=dict(value=...), "
             "smooth_width=dict(value=...))")
        .def("f_c", py::vectorize(&CutoffFunction_t::f_c))
        .def("df_c", py::vectorize(&CutoffFunction_t::df_c), "df_c(r) / dr");

    using CutoffFunction_rs_t =
        internal::CutoffFunction<internal::CutoffFunctionType::RadialScaling>;

    py::class_<CutoffFunction_rs_t> fc_rs(mod,
                                          "CutoffFunctionShiftedRadialScaling");
    fc_rs
        .def(py::init([](const py::dict & hyper) {
               // convert to json
               json hypers = hyper;
               return std::make_unique<CutoffFunction_rs_t>(hypers);
             }),
             R"(hyper should be a
           dict(cutoff=dict(value=XXX),
                smooth_width=dict(value=XXX),
                rate=dict(value=XXX),
                exponent=dict(value=XXX),
                scale=dict(value=XXX)))")
        .def("f_c", py::vectorize(&CutoffFunction_rs_t::f_c))
        .def("df_c", py::vectorize(&CutoffFunction_rs_t::df_c), "df_c(r) / dr")
        .def_property_readonly("cutoff",
            [](CutoffFunction_rs_t & fc) {
              return fc.cutoff;
            });
  }

  void add_math(py::module & mod) {
    bind_sph(mod);
    bind_ri(mod);
    bind_fc(mod);
  }
}  // namespace rascal
