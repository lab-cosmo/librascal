/**
 * @file   bind_include.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   31 Oct 2018
 *
 * @brief  File to centralize includes and function declaration
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

#ifndef BINDINGS_BIND_INCLUDE_HH_
#define BINDINGS_BIND_INCLUDE_HH_

#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/utils/json_io.hh"
#include "rascal/utils/utils.hh"

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
// for the hasattr function to test the module namespace
#include <pybind11/pytypes.h>

#include <Eigen/Dense>

#include <map>
#include <memory>
#include <vector>

/*
 * Prevent vector of atomic structures from being copied into a Python list,
 * since we already have the AtomsList object.  See also
 * https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html
 */
PYBIND11_MAKE_OPAQUE(std::vector<rascal::AtomicStructure<3>>);

namespace py = pybind11;

//! Simplistic but robust implicit conversion of py::dict to/from nlohmann::json
namespace nlohmann {
  template <>
  struct adl_serializer<py::dict> {
    static void to_json(json & j, const py::dict & dic) {
      py::module py_json = py::module::import("json");
      j = json::parse(
          static_cast<std::string>(py::str(py_json.attr("dumps")(dic))));
    }
    static void from_json(const json & j, py::dict & dic) {
      py::module py_json = py::module::import("json");
      dic = py_json.attr("loads")(j.dump());
    }
  };
}  // namespace nlohmann

namespace rascal {
  namespace internal {
    /**
     * Expose to python the serialization of Object as a python dictionary.
     *
     * @tparam Object is expected to be nlohmann::json (de)serializable
     *
     * A copy and a json (de)serialization are necessary to make sure that if
     * the resulting dictionary is written in json, then it will be directly
     * convertible to the original object in C++ and vice-versa.
     */
    template <class Object, class... Bases>
    void bind_dict_representation(py::class_<Object, Bases...> & obj) {
      // serialization to a python dictionary
      obj.def("to_dict", [](const Object & self) {
        json j;
        j = self;  // implicit conversion to nlohmann::json
        return j.template get<py::dict>();
      });
      // construction from a python dictionary
      obj.def_static("from_dict", [](const py::dict & d) {
        json j;
        j = d;  // implicit conversion to nlohmann::json
        return std::make_unique<Object>(j.template get<Object>());
      });
      // string representation
      obj.def("__str__", [](const Object & self) {
        json j = self;  // implicit conversion to nlohmann::json
        std::string str = j.dump(2);
        std::string representation_name{internal::type_name<Object>()};
        std::string sep{" | Parameters: "};
        std::string prefix{"Class: "};
        return prefix + representation_name + sep + str;
      });
    }

    /**
     * Transforms the template type to a string for the python bindings.
     * There are submodules in the python bindings with the class
     * title so to avoid redundancy they are removed from the
     * typename.
     * @tparam T type that should be stringified
     * @returns std::string name of the type
     */
    template <typename T>
    std::string GetBindingTypeName() {
      static std::map<std::string, std::string> replacement_map = {
          {"StructureManager", ""}, {"Adaptor", ""}, {"Calculator", ""}};

      std::string name = type_name<T>();
      for (const auto & map : replacement_map) {
        replace(name, map.first, map.second);
      }

      return name;
    }
  }  // namespace internal
}  // namespace rascal

#endif  // BINDINGS_BIND_INCLUDE_HH_
