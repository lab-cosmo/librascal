/**
 * file   math_interface.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief defines an interface to other math library like cephes to separate the
 *        namespaces
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "math_interface.hh"
#include "math/cephes/mconf.h"

namespace rascal {
  namespace math {

    double hyp2f1(const double & a, const double & b, const double & c,
                  const double & x) {
      return cephes::hyp2f1(a, b, c, x);
    }

    double hyp1f1(const double & a, const double & b, const double & x) {
      return cephes::hyperg(a, b, x);
    }

    double pow_u(double x, size_t n) {
      double value = 1.0;

      /* repeated squaring method
       * returns 0.0^0 = 1.0, so continuous in x
       * (from GSL)
       */
      do {
        if (n & 1)
          value *= x; /* for n odd */
        n >>= 1;
        x *= x;
      } while (n);

      return value;
    }

    double pow_i(double x, int n) {
      size_t un;

      if (n < 0) {
        x = 1.0 / x;
        un = static_cast<size_t>(-n);
      } else {
        un = static_cast<size_t>(n);
      }

      return pow_u(x, un);
    }

    //    double pow_cephes(const double & x, const int & n) { return
    //    cephes::powi(x, n); }
    /*
        double pow(const double & x, const std::size_t & n) {
          return cephes::powi(x, static_cast<int>(n));
        }
    */
  }  // namespace math
}  // namespace rascal
