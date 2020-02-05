/**
 * @file   rascal/utils/units.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief  unit standardisation
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

#ifndef SRC_RASCAL_UTILS_UNITS_HH_
#define SRC_RASCAL_UTILS_UNITS_HH_

#include <map>
#include <string>

namespace rascal {

  namespace units {

    class UnitStyle {
     public:
      //! Default constructor
      UnitStyle() = delete;

      //! full constructor
      UnitStyle(const std::string & mass, const std::string & distance,
                const std::string & time, const std::string & energy,
                const std::string & velocity, const std::string & force,
                const std::string & torque, const std::string & temperature,
                const std::string & pressure,
                const std::string & dynamic_viscosity,
                const std::string & charge, const std::string & dipole,
                const std::string & electric_field,
                const std::string & density);

      //! copy constructor
      UnitStyle(const UnitStyle & other) = default;

      //! Move constructor
      UnitStyle(UnitStyle && other) = default;

      //! Destructor
      virtual ~UnitStyle() = default;

      //! Copy assignment operator
      UnitStyle & operator=(const UnitStyle & other) = delete;

      //! Move assignment operator
      UnitStyle & operator=(UnitStyle && other) = delete;

      /**
       * returns a string representation of the given unit with a optional
       * fractional exponent. E.g., in the metal system `mass(2)` returns
       * '(grams/mole)^(2)', `temperature(-1, 2)` returns '(K)^(-1/2)', or a
       * `simple distance()` returns 'Å'.
       */
      const std::string mass(int numerator = 1, int denominator = 1) const;
      const std::string distance(int numerator = 1, int denominator = 1) const;
      const std::string time(int numerator = 1, int denominator = 1) const;
      const std::string energy(int numerator = 1, int denominator = 1) const;
      const std::string velocity(int numerator = 1, int denominator = 1) const;
      const std::string force(int numerator = 1, int denominator = 1) const;
      const std::string torque(int numerator = 1, int denominator = 1) const;
      const std::string temperature(int numerator = 1,
                                    int denominator = 1) const;
      const std::string pressure(int numerator = 1, int denominator = 1) const;
      const std::string dynamic_viscosity(int numerator = 1,
                                          int denominator = 1) const;
      const std::string charge(int numerator = 1, int denominator = 1) const;
      const std::string dipole(int numerator = 1, int denominator = 1) const;
      const std::string electric_field(int numerator = 1,
                                       int denominator = 1) const;
      const std::string density(int numerator = 1, int denominator = 1) const;

      const std::string none(int numerator = 1, int denominator = 1) const;

     protected:
      const std::string format(const std::string & symbol, int numerator,
                               int denominator) const;
      const std::string _mass;
      const std::string _distance;
      const std::string _time;
      const std::string _energy;
      const std::string _velocity;
      const std::string _force;
      const std::string _torque;
      const std::string _temperature;
      const std::string _pressure;
      const std::string _dynamic_viscosity;
      const std::string _charge;
      const std::string _dipole;
      const std::string _electric_field;
      const std::string _density;
      const std::string _none{"-"};
    };

    extern const UnitStyle metal;
    extern const UnitStyle electron;

    extern const std::map<std::string, int> species_numbers;

  }  // namespace units

}  // namespace rascal

#endif  // SRC_RASCAL_UTILS_UNITS_HH_
