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

#ifndef SRC_MATH_MATH_INTERFACE_HH_
#define SRC_MATH_MATH_INTERFACE_HH_

#include <cmath>

/*
 * Cephes polutes the namespace with extern definitions so we don't link to it
 * directly. math_interface.cc includes math/cephes/mconf.h which in turns
 * includes math/cephes/cephes.h (the main header with all the functions). It
 * also redefines the functions needed, e.g. hyp2f1, which is exposed in
 * math_interface.hh. So when calling rascal::math::hyp2f1 you call a function
 * from librascal which uses a function from libcephes. That is how the
 * namespace is not polluted because only math_interface.cc sees it and it is a
 * .cc
 */

namespace rascal {
  namespace math {

    double hyp2f1(const double & a, const double & b, const double & c,
                  const double & x);
    double hyp1f1(const double & a, const double & b, const double & x);

    // overload the pow function from cmath to optimize for integer power
    double pow(const double & x, const int & n);
    double pow(const double & x, const std::size_t & n);

    inline double pow(const double & x, const double & n) {
      return std::pow(x, n);
    }
  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_MATH_INTERFACE_HH_
