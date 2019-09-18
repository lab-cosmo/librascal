/**
 * @file   models.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   16 Jun 2019
 *
 * @brief  File for binding the models
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "bind_py_models.hh"

namespace rascal {

  decltype(auto) add_kernel(py::module & mod, py::module & /*m_internal*/) {
    py::class_<Kernel> kernel(mod, "Kernel");
    // use custom constructor to pass json formated string as initializer
    // an alternative would be to convert python dict to json internally
    // but needs some work on in the pybind machinery
    kernel.def(py::init([](std::string & hyper_str) {
      // convert to json
      json hypers = json::parse(hyper_str);
      return std::make_unique<Kernel>(hypers);
    }));

    return kernel;
  }

  template <internal::KernelType Type, class Calculator,
            class StructureManagers, class CalculatorBind>
  void bind_kernel_compute_function(CalculatorBind & kernel) {
    kernel.def("compute",
               &Kernel::template compute<Calculator, StructureManagers>,
               py::call_guard<py::gil_scoped_release>());
  }

  /**
   * Function to bind the representation managers to python
   *
   * @params mod pybind11 representation of the python module the represenation
   *             managers will be included to
   * @params m_internal pybind11 representation of the python module that are
   *                  needed but not useful to use on the python side
   *
   */
  void add_kernels(py::module & mod, py::module & m_internal) {
    // Defines a particular structure manager type
    using ManagerCollection_t =
        ManagerCollection<StructureManagerCenters, AdaptorNeighbourList,
                          AdaptorStrict>;

    // Defines the representation manager type for the particular structure
    // manager
    using Calc1_t = CalculatorSphericalInvariants;
    // Bind the interface of this representation manager
    auto kernel = add_kernel(mod, m_internal);
    bind_kernel_compute_function<internal::KernelType::Cosine, Calc1_t,
                                 ManagerCollection_t>(kernel);
  }
}  // namespace rascal
