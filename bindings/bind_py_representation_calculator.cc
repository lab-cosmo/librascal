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

#include "bind_py_representation_calculator.hh"

namespace rascal {

  template <class Calculator>
  using PyCalculator = py::class_<Calculator, CalculatorBase>;

  template <typename Calculator>
  auto add_representation_calculator(py::module & mod,
                                     py::module & /*m_unused*/) {
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

    return representation;
  }

  template <class Calculator, typename Manager,
            template <class> class... Adaptor>
  void bind_compute_function(PyCalculator<Calculator> & representation) {
    using TypeHolder_t = StructureManagerTypeHolder<Manager, Adaptor...>;
    using Manager_t = typename TypeHolder_t::type;
    representation.def(
        "compute", &Calculator::template compute<std::shared_ptr<Manager_t>>,
        py::call_guard<py::gil_scoped_release>());

    representation.def(
        "compute",
        &Calculator::template compute<ManagerCollection<Manager, Adaptor...>>,
        py::call_guard<py::gil_scoped_release>());
  }

  namespace py_internal {
    template <typename SM, typename AdaptorTypeHolder_>
    struct bind_compute_function_helper;

    template <typename SM, template <class> class... Ti>
    struct bind_compute_function_helper<SM, AdaptorTypeHolder<Ti...>> {
      template <class Calculator>
      static void apply(PyCalculator<Calculator> & representation) {
        bind_compute_function<Calculator, SM, Ti...>(representation);
      }
    };

    template <typename StructureManagerTypeHolder_>
    struct bind_compute_function_util;

    template <typename... T>
    struct bind_compute_function_util<std::tuple<T...>> {
      template <class Calculator>
      static void apply(PyCalculator<Calculator> & representation) {
        bind_compute_function_helper<T...>::apply(representation);
      }
    };
  }  // namespace py_internal

  template <typename StructureManagerTypeHolder_, class Calculator>
  void bind_compute_function_helper(PyCalculator<Calculator> & representation) {
    py_internal::bind_compute_function_util<StructureManagerTypeHolder_>::apply(
        representation);
  }

  /**
   * Function to bind the representation managers to python
   *
   * @param mod pybind11 representation of the python modules the representation
   *        managers will be included to
   * @param m_internal pybind11 representation of the python modules that are
   *        needed but not useful to use on the python side
   */
  void add_representation_calculators(py::module & mod,
                                      py::module & m_internal) {
    auto base = py::class_<CalculatorBase>(m_internal, "CalculatorBase");
    base.def_readwrite("name", &CalculatorBase::name);
    base.def_readonly("default_prefix", &CalculatorBase::default_prefix);
    /*-------------------- rep-bind-start --------------------*/
    // Defines a particular structure manager type

    using TypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;
    using ManagerList_t = typename TypeHolder_t::type_list;
    // using Manager_t = typename TypeHolder_t::type;
    // StructureManagerCenters,AdaptorNeighbourList, AdaptorStrict
    // Defines the representation manager type for the particular structure
    // manager
    using Calc1_t = CalculatorSortedCoulomb;
    // Bind the interface of this representation manager
    auto rep_sorted_coulomb =
        add_representation_calculator<Calc1_t>(mod, m_internal);
    bind_compute_function_helper<ManagerList_t>(rep_sorted_coulomb);

    /*-------------------- rep-bind-end --------------------*/
    using ManagerList_1_t = typename StructureManagerTypeHolder<
        StructureManagerCenters, AdaptorNeighbourList,
        AdaptorCenterContribution, AdaptorStrict>::type_list;
    using Calc2_t = CalculatorSphericalExpansion;
    auto rep_spherical_expansion =
        add_representation_calculator<Calc2_t>(mod, m_internal);
    bind_compute_function_helper<ManagerList_1_t>(rep_spherical_expansion);

    using Calc3_t = CalculatorSphericalInvariants;
    auto rep_soap = add_representation_calculator<Calc3_t>(mod, m_internal);
    bind_compute_function_helper<ManagerList_1_t>(rep_soap);

    using Calc4_t = CalculatorSphericalCovariants;
    auto rep_lambda_soap =
        add_representation_calculator<Calc4_t>(mod, m_internal);
    bind_compute_function_helper<ManagerList_1_t>(rep_lambda_soap);
  }

}  // namespace rascal
