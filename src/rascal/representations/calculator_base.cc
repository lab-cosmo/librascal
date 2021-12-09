/**
 * @file   rascal/representations/calculator_base.cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  base class for representation managers
 *
 * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "rascal/representations/calculator_base.hh"

namespace rascal {

  std::string CalculatorBase::get_options_string() {
    std::string out{};
    for (const auto & item : this->options) {
      out += item.second + std::string(" ");
    }
    if (not out.empty()) {
      out.pop_back();
    }
    return out;
  }
  /* ---------------------------------------------------------------------- */
  std::string CalculatorBase::get_hypers_string() {
    return this->hypers.dump(2);
  }
  /* ---------------------------------------------------------------------- */
  void CalculatorBase::check_hyperparameters(
      const CalculatorBase::ReferenceHypers_t & reference_items,
      const CalculatorBase::Hypers_t & hypers) {
    for (const auto & reference_item : reference_items) {
      auto && key{reference_item.first};
      auto && val{reference_item.second};

      // check if the registered key belongs to the input dictionary

      if (hypers.find(key) != hypers.end()) {
        // if there are few possible values for the value of reference_item
        // check if at least one is in hypers
        if (not val.empty()) {
          bool is_valid_value{false};
          for (const auto & value : val) {
            if (value == hypers[key].get<std::string>()) {
              is_valid_value = true;
            }
          }
          // if the value does not appear in the list of possible values
          if (not is_valid_value) {
            auto error_message{std::string("Option ") + key +
                               std::string(" with value ") +
                               hypers[key].get<std::string>() +
                               std::string(" is not implemented")};
            throw std::invalid_argument(error_message.c_str());
          }
        }
      } else {
        auto error_message{std::string("Parameter '") + key +
                           std::string("' is missing from the inputs.")};
        throw std::invalid_argument(error_message.c_str());
      }
    }
  }

}  // namespace rascal
