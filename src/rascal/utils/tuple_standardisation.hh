/**
 * file   tuple_standardisation.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief helper class to create a common type representing std::arrays of
 * varying length for use as keys
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

#ifndef SRC_RASCAL_UTILS_TUPLE_STANDARDISATION_HH_
#define SRC_RASCAL_UTILS_TUPLE_STANDARDISATION_HH_

#include <array>
#include <limits>

namespace rascal {

  template <typename T, size_t MaxOrder>
  class TupleStandardisation {
    static_assert(std::is_arithmetic<T>::value,
                  "must be an integer or float type");

   public:
    //! Default constructor
    TupleStandardisation() = delete;

    template <size_t Order>
    explicit TupleStandardisation(const std::array<T, Order> & tup)
        : order{Order} {
      for (size_t i{0}; i < Order; ++i) {
        this->key[i] = std::numeric_limits<T>::min();
      }
      for (size_t i{0}; i < MaxOrder - Order; ++i) {
        this->key[i + Order] = tup[i];
      }
    }

    //! Copy constructor
    TupleStandardisation(const TupleStandardisation & other) = default;

    //! Move constructor
    TupleStandardisation(TupleStandardisation && other) = default;

    //! Destructor
    virtual ~TupleStandardisation() = default;

    //! Copy assignment operator
    TupleStandardisation &
    operator=(const TupleStandardisation & other) = default;

    //! Move assignment operator
    TupleStandardisation & operator=(TupleStandardisation && other) = default;

    bool operator<(const TupleStandardisation & other) const {
      return this->key < other.key;
    }

    template <size_t Order>
    operator std::array<T, Order>() const {
      std::array<T, Order> ret_val{};
      for (size_t i{0}; i < Order; ++i) {
        ret_val[i] = this->key[i + MaxOrder - Order];
      }
      return ret_val;
    }

    T & at(size_t index) { return this->key.at(index); }
    const T & at(size_t index) const { return this->key.at(index); }
    const size_t & get_order() const { return this->order; }

   protected:
    const std::array<T, MaxOrder> key{};
    const size_t order;
  };

}  // namespace rascal

#endif  // SRC_RASCAL_UTILS_TUPLE_STANDARDISATION_HH_
