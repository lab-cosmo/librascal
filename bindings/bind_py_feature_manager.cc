/**
 * @file   bind_py_feature_manager.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   30 Oct 2018
 *
 * @brief  File for binding the Feature Managers
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

using namespace rascal; // NOLINT

/**
 * Bind a feature manager
 */
template<template<class> class FeatureManager_t, typename T>
decltype(auto) bind_feature_manager(py::module & mod, py::module & ) {
  using Feature = FeatureManager_t<T>;

  std::string feature_name =
        internal::GetBindingTypeName<Feature>();

  py::class_<Feature, FeatureManagerBase<T> >
             feature(mod, feature_name.c_str());
  // use custom constructor to pass json formated string as initializer
  // an alternative would be to convert python dict to json internally
  // but needs some workon in the pybind machinery
  feature.def(py::init([](int & n_feature, std::string& hyper_str) {
        // convert to json
        json hypers = json::parse(hyper_str);
        return std::make_unique<Feature>(n_feature, hypers);
            }));
  feature.def("reserve", &Feature::reserve);
  feature.def("append",
        (void (Feature::*)(RepresentationManagerBase&)) &Feature::push_back);
  feature.def_property_readonly("size", &Feature::size,
                              py::return_value_policy::copy);
  feature.def_property_readonly("shape", &Feature::shape,
                              py::return_value_policy::copy);
  feature.def("get_feature_matrix", &Feature::get_feature_matrix,
        py::return_value_policy::reference_internal, py::keep_alive<1, 0>());

  return feature;
}

//! Feature aggregator python binding
void add_feature_managers(py::module & mod, py::module & m_garbage) {
  using FeatureBase_0 = FeatureManagerBase<double>;
  std::string featurebase_0_name =
        internal::GetBindingTypeName<FeatureBase_0>();
  py::class_<FeatureManagerBase<double>>(m_garbage,
                                          featurebase_0_name.c_str());
  auto feature_double =
         bind_feature_manager<FeatureManagerDense, double>(mod, m_garbage);

  using FeatureBase_1 = FeatureManagerBase<float>;
  std::string featurebase_1_name =
        internal::GetBindingTypeName<FeatureBase_1>();
  py::class_<FeatureManagerBase<float>>(m_garbage,
                                        featurebase_1_name.c_str());
  auto feature_float =
      bind_feature_manager<FeatureManagerDense, float>(mod, m_garbage);
}
