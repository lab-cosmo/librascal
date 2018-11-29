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
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "math_interface.hh"
#include "math/cephes/mconf.h"


namespace rascal {
  namespace math {

    double hyp2f1(const double & a, const double & b, const double & c,
                  const double & x );
      return cephes::hyp2f1(a, b, c, x);
    }

    double hyp1f1(const double& a, const double& b, const double& x) {
      return cephes::hyperg(a, b, x);
    }

  } // math
} // rascal
