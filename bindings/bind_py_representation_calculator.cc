/**
 * @file   bind_py_representation_calculator.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   30 Oct 2018
 *
 * @brief  File for binding the Representation calculator
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "bind_include.hh"

using namespace rascal;  // NOLINT

template <typename Calculator>
decltype(auto) add_representation_calculator(py::module & mod, py::module &) {

  std::string representation_name =
      internal::GetBindingTypeName<Calculator>();

  py::class_<Calculator, CalculatorBase> representation(
      mod, representation_name.c_str());
  // use custom constructor to pass json formated string as initializer
  // an alternative would be to convert python dict to json internally
  // but needs some work on in the pybind machinery
  representation.def(py::init([](std::string & hyper_str) {
    // convert to json
    json hypers = json::parse(hyper_str);
    return std::make_unique<Calculator>(hypers);
  }));

  representation.def("get_name",
                     &Calculator::get_name,
                     py::call_guard<py::gil_scoped_release>());

  return representation;
}

template <class Calculator, class StructureManager, class CalculatorBind>
void bind_compute_function(CalculatorBind& representation) {
  representation.def("compute",
                     &Calculator::template compute<std::shared_ptr<StructureManager>>,
                     py::call_guard<py::gil_scoped_release>());
}

/**
 * Function to bind the representation managers to python
 *
 * @params mod pybind11 representation of the python module the represenation
 *             managers will be included to
 * @params m_throwaway pybind11 representation of the python module that are
 *                  needed but not useful to use on the python side
 *
 */
void add_representation_calculators(py::module & mod, py::module & m_throwaway) {
  py::class_<CalculatorBase>(m_throwaway, "CalculatorBase");
  /*-------------------- rep-bind-start --------------------*/
  // Defines a particular structure manager type
  using Manager_t =
      AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>;
  // Defines the representation manager type for the particular structure
  // manager
  using Calc1_t = CalculatorSortedCoulomb;
  // Bind the interface of this representation manager
  auto rep_sorted_coulomb =
      add_representation_calculator<Calc1_t>(mod, m_throwaway);
  bind_compute_function<Calc1_t, Manager_t>(rep_sorted_coulomb);
  /*-------------------- rep-bind-end --------------------*/
  using Calc2_t = CalculatorSphericalExpansion;
  auto rep_spherical_expansion =
      add_representation_calculator<Calc2_t>(mod, m_throwaway);
  bind_compute_function<Calc2_t, Manager_t>(rep_spherical_expansion);

  using Calc3_t = CalculatorSphericalInvariants;
  auto rep_soap =
      add_representation_calculator<Calc3_t>(mod, m_throwaway);
  bind_compute_function<Calc3_t, Manager_t>(rep_soap);
}
