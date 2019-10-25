/**
 * @file   key_standardisation.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief helper class to create a common type representing std::arrays of
 * varying length for use as keys in maps
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

#ifndef SRC_UTILS_KEY_STANDARDISATION_HH_
#define SRC_UTILS_KEY_STANDARDISATION_HH_

#include <array>
#include <limits>

namespace rascal {

  template <typename T, size_t MaxOrder>
  class KeyStandardisation {
    static_assert(std::is_arithmetic<T>::value,
                  "must be an integer or float type");

   public:
    //! Default constructor
    KeyStandardisation() = delete;

    template <size_t Order>
    explicit KeyStandardisation(const std::array<T, Order> & tup)
        : key{KeyStandardisation::fill(tup)}, order{Order} {}

    //! Copy constructor
    KeyStandardisation(const KeyStandardisation & other) = default;

    //! Move constructor
    KeyStandardisation(KeyStandardisation && other) = default;

    //! Destructor
    virtual ~KeyStandardisation() = default;

    //! Copy assignment operator
    KeyStandardisation & operator=(const KeyStandardisation & other) = default;

    //! Move assignment operator
    KeyStandardisation & operator=(KeyStandardisation && other) = default;

    bool operator<(const KeyStandardisation & other) const {
      return this->key < other.key;
    }

    template <size_t Order>
    operator const std::array<T, Order> &() const {
      static_assert(Order <= MaxOrder,
                    "Casting a larger than allowed key size");
      return reinterpret_cast<std::array<T, Order> &>(*this);
    }

    const T & operator[](const size_t index) { return this->key[index]; }

    const T & operator[](const size_t index) const { return this->key[index]; }

    size_t get_order() const { return this->order; }

    const T & back() const { return key.back(); }

   protected:
    const std::array<T, MaxOrder> key;
    const size_t order;

    //! function for filling the key of the standardised tuple, which in
    //! practise is an array
    template <size_t Order>
    static std::array<T, MaxOrder> fill(const std::array<T, Order> & tup) {
      static_assert(Order <= MaxOrder,
                    "Impossible key type, Order has maximum value of MaxOrder");
      std::array<T, MaxOrder> key;
      for (size_t i{Order}; i < MaxOrder; ++i) {
        key[i] = std::numeric_limits<T>::min();
      }
      for (size_t i{0}; i < Order; ++i) {
        key[i] = tup[i];
      }
      return key;
    }
  };
}  // namespace rascal

#endif  // SRC_UTILS_KEY_STANDARDISATION_HH_
