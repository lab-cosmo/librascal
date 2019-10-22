/**
 * @file   json_io.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief  implementations of json manipulation functions
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

#include "json_io.hh"

namespace rascal {
  namespace json_io {

    json load(const std::string & filename) {
      json data;
      auto extension{internal::get_filename_extension(filename)};
      if (extension == "json") {
        data = load_txt(filename);
      } else if (extension == "ubjson") {
        data = load_bin(filename);
      } else {
        throw std::runtime_error(std::string("Don't know the extension of ") +
                                 filename);
      }
      return data;
    }

    json load_txt(const std::string & filename) {
      json j;
      std::ifstream reader(filename);
      if (not reader.is_open()) {
        throw std::runtime_error(std::string("Could not open the file: ") +
                                 filename);
      }
      reader >> j;
      reader.close();
      return j;
    }

    json load_bin(const std::string & filename) {
      return json::from_ubjson(internal::read_binary_file(filename));
    }

    /* ---------------------------------------------------------------------- */
    void to_json(json & j, AtomicJsonData & s) {
      j = json{{"cell", s.cell},
               {"atom_types", s.type},
               {"pbc", s.pbc},
               {"positions", s.position}};
    }

    /* ---------------------------------------------------------------------- */
    void from_json(const json & j, AtomicJsonData & s) {
      if (j.count("atom_types") == 1) {
        s.type = j.at("atom_types").get<std::vector<int>>();
      } else if (j.count("numbers") == 1) {
        s.type = j.at("numbers").get<std::vector<int>>();
      } else {
        throw std::runtime_error(
            R"(AtomicJsonData needs atom_types or numbers keyword)");
      }
      s.cell = j.at("cell").get<std::vector<std::vector<double>>>();
      s.pbc = j.at("pbc").get<std::vector<int>>();
      s.position = j.at("positions").get<std::vector<std::vector<double>>>();
    }

    /* ---------------------------------------------------------------------- */
    double check_units(const std::string & expected_unit,
                       const json & parameter) {
      if (not(expected_unit == parameter.at("unit").get<std::string>())) {
        std::stringstream error{};
        error << "unit '" << parameter.at("unit").get<std::string>()
              << "' differs from the expected unit '" << expected_unit << "'.";
        throw std::runtime_error(error.str());
      }
      return parameter.at("value").get<double>();
    }

  }  // namespace json_io

}  // namespace rascal
