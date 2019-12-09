/**
 * @file   rascal/units.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief  implementation for unit standardisation
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

#include "rascal/units.hh"

#include <sstream>

namespace rascal {

  namespace units {

    /* ---------------------------------------------------------------------- */
    UnitStyle::UnitStyle(const std::string & mass, const std::string & distance,
                         const std::string & time, const std::string & energy,
                         const std::string & velocity,
                         const std::string & force, const std::string & torque,
                         const std::string & temperature,
                         const std::string & pressure,
                         const std::string & dynamic_viscosity,
                         const std::string & charge, const std::string & dipole,
                         const std::string & electric_field,
                         const std::string & density)
        : _mass{mass}, _distance{distance}, _time{time}, _energy{energy},
          _velocity{velocity}, _force{force}, _torque{torque},
          _temperature{temperature}, _pressure{pressure},
          _dynamic_viscosity{dynamic_viscosity}, _charge{charge},
          _dipole{dipole}, _electric_field{electric_field}, _density{density} {}

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::mass(int numerator, int denominator) const {
      return this->format(this->_mass, numerator, denominator);
    }
    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::distance(int numerator,
                                          int denominator) const {
      return this->format(this->_distance, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::time(int numerator, int denominator) const {
      return this->format(this->_time, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::energy(int numerator, int denominator) const {
      return this->format(this->_energy, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::velocity(int numerator,
                                          int denominator) const {
      return this->format(this->_velocity, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::force(int numerator, int denominator) const {
      return this->format(this->_force, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::torque(int numerator, int denominator) const {
      return this->format(this->_torque, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::temperature(int numerator,
                                             int denominator) const {
      return this->format(this->_temperature, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::pressure(int numerator,
                                          int denominator) const {
      return this->format(this->_pressure, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::dynamic_viscosity(int numerator,
                                                   int denominator) const {
      return this->format(this->_dynamic_viscosity, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::charge(int numerator, int denominator) const {
      return this->format(this->_charge, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::dipole(int numerator, int denominator) const {
      return this->format(this->_dipole, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::electric_field(int numerator,
                                                int denominator) const {
      return this->format(this->_electric_field, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::density(int numerator, int denominator) const {
      return this->format(this->_density, numerator, denominator);
    }

    /* ---------------------------------------------------------------------- */
    const std::string UnitStyle::format(const std::string & symbol,
                                        int numerator, int denominator) const {
      if ((numerator == 1) && (denominator == 1)) {
        return symbol;
      } else if (denominator == 1) {
        std::stringstream out_stream{};
        out_stream << "(" << symbol << ")^(" << numerator << ")";
        return out_stream.str();
      } else {
        std::stringstream out_stream{};
        out_stream << "(" << symbol << ")^(" << numerator << "/" << denominator
                   << ")";
        return out_stream.str();
      }
    }

    /* ---------------------------------------------------------------------- */

    static const std::string undefined{"undefined"};
    const  UnitStyle metal{"(g/mol)", "Å",     "ps",    "eV",          "(Å/ps)",
                          "(eV/Å)",  "eV",    "K",     "bar",         "P",
                          "e",       "(e*Å)", "(V/Å)", "((g/cm)^(3))"};
    const UnitStyle electron{
        "u", "rBohr", "fs",      "Ha", "(rBohr/atu)", "(Ha/rBohr)", undefined,
        "K", "Pa",    undefined, "e",  "D",           "(V/cm)",     undefined};

    /* ---------------------------------------------------------------------- */
    namespace internal {
      std::map<std::string, int> make_species() {
        std::map<std::string, int> ret_val{};
        ret_val["Mg"] = 12;
        ret_val["Si"] = 14;
        return ret_val;
      }
    }  // internal

    /* ---------------------------------------------------------------------- */
    const std::map<std::string, int> species_numbers{internal::make_species()};

  }  // namespace units

}  // namespace rascal
