/**
 * file   enum_map.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   18 March 2019
 *
 * @brief
 *
 * Copyright © 2018 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_ENUM_MAP_HH_
#define SRC_ENUM_MAP_HH_

#include <array>
#include <exception>
#include <string>

namespace rascal {
  namespace internal {

    // I like exceptions, you might feel differently
    class UnknownValueException : public std::runtime_error {
    public:
      UnknownValueException(const std::string& name):std::runtime_error("Unknown value: " + name) {};
      UnknownValueException(int value):std::runtime_error("Unknown name for enum value: " + std::to_string(value)) {};
    };

    template<class T>
    struct NameValuePair {
      using value_type = T;
      const T value;
      const char* const name;
    };

    // Templated helper functions.
    // Mapping is some type of standard container that supports find_if()
    // V is the type of the enum whose value we wish to look up
    template<class Mapping, class V>
    std::string getNameForValue(Mapping a, V value) {
      auto pos = std::find_if(std::begin(a), std::end(a), [&value](const typename Mapping::value_type& t){
          return (t.value == value);
      });
      if (pos != std::end(a)) {
          return pos->name;
      }

      throw UnknownValueException(static_cast<int>(value));
      // or return some default value here
      // return Mapping::value_type::value_type();
    }

    template<class Mapping>
    typename Mapping::value_type::value_type getValueForName(Mapping a, const std::string& name) {
      auto pos = std::find_if(std::begin(a), std::end(a), [&name](const typename Mapping::value_type& t){
          return (t.name == name);
      });
      if (pos != std::end(a)) {
          return pos->value;
      }

      throw UnknownValueException(name);
      // or return an empty string, whatever works for you
    }


  } // namespace internal
} // namespace rascal

#endif  // SRC_ENUM_MAP_HH_
