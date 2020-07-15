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

#include "rascal/math/spherical_harmonics.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"

namespace rascal {

  void bind_sph(py::module & mod) {
    py::class_<SphericalHarmonics> sph(
        mod, "SphericalHarmonics");
    sph.def("compute",
              py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(&SphericalHarmonics::calc)
    );
    sph.def("harmonics", &SphericalHarmonics::get_harmonics);
    sph.def("harmonics_derivatives",
            &SphericalHarmonics::get_harmonics_derivatives);
  }

  void bind_ri(py::module & mod) {
    py::class_<RadialContribution<RadialBasisType::DVR>> ri_dvr(
        mod, "RadialContributionDVR");
    ri_dvr.def(
      "",
    );
    ri_dvr.def(py::init([](const py::dict & hyper) {
      // convert to json
      json hypers = hyper;
      return std::make_unique<RadialContribution<RadialBasisType::DVR>>(hypers);
    }));

    ri_dvr.def("compute_center_contribution",
              py::overload_cast<double>(&RadialContribution::template <RadialBasisType::DVR>::compute_center_contribution)
    );

    ri_dvr.def("compute_neighbour_contribution",
              py::overload_cast<const double, const double>(&RadialContribution::template <RadialBasisType::DVR>::compute_neighbour_contribution)
    );

    ri_dvr.def("compute_neighbour_derivative",
              py::overload_cast<>(&RadialContribution::template <RadialBasisType::DVR>::compute_neighbour_derivative)
    );


  }

  void math_binding(py::module & mod) {
    py::module m_math = mod.def_submodule("math");
    m_math.doc() = "Math Routines";
    bind_sph(m_math);
    bind_ri(m_math);
  }
}  // namespace rascal
