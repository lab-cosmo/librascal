/**
 * file   json_io.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   18 Jun 2018
 *
 * @brief JSON interface from nlohmanns header class
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef JSON_IO_H
#define JSON_IO_H

// An external header-library/header-class, which makes it easy to use
// the JSON as a first class data type. See
// https://github.com/nlohmann/json for documentation.
#include "json.hpp"

// For convenience
using json = nlohmann::json;

// All functions and classes are in the namespace <code>rascal</code>,
// which ensures that they don't clash with other libraries one might
// use in conjunction.
namespace rascal {
  namespace json_io {
    // To read from a JSON file and deserialize the content, the used
    // class needs a <code>struct</code> with standard data types.
    struct AtomicStructure {
      /**
         \param cell is a vector a vector of vectors which holds the cell unit
         vectors.
         \param type a vector of integers which holds the atomic type
         (coordination number).
         \param pbc is a 0/1 vector which says, where periodic boundary
         conditions are applied.
         \param position is a vector of vectors which holds the atomic
         positions.
      */
      std::vector<std::vector<double>> cell{};
      std::vector<int> type{};
      std::vector<int> pbc{};
      std::vector<std::vector<double>> position{};
    };

    // This function is used to convert to the JSON format with the
    // given keywords. It is an overload of the function defined in the header
    // class json.hpp.
    // Inline needed, otherwise it is a multiple definition
    inline void to_json(json & j, AtomicStructure& s) {
      j = json{
        {"cell", s.cell},
        {"numbers", s.type},
        {"pbc", s.pbc},
        {"positions", s.position}
      };
    }

    // This function is used to read from the JSON file and convert
    // the data into standard types. It is an overload of the function defined
    // in json.hpp class header.
    inline void from_json(const json& j, AtomicStructure& s) {
      s.cell = j.at("cell").get<std::vector<std::vector<double>>>();
      s.type = j.at("numbers").get<std::vector<int>>();
      s.pbc = j.at("pbc").get<std::vector<int>>();
      s.position = j.at("positions").get<std::vector<std::vector<double>>>();
    }
  }
} // rascal


#endif /* JSON_IO_H */
