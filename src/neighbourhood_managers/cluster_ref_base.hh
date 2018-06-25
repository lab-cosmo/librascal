/**
 * file   ClusterRefBase.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   21 Jun 2018
 *
 * @brief  a base class for getting access to clusters
 *
 * Copyright © 2018 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef CLUSTERREFBASE_H
#define CLUSTERREFBASE_H

#include <array>

namespace rascal {

  template<int Level>
  class ClusterRefBase
  {
  public:
    //! Default constructor
    ClusterRefBase() = delete;

    ClusterRefBase(std::array<int, Level> indices):
    indices{indices} {}

    //! Copy constructor
    ClusterRefBase(const ClusterRefBase &other) = default; //

    //! Move constructor
    ClusterRefBase(ClusterRefBase &&other) = default;

    //! Destructor
    virtual ~ClusterRefBase() = default;

    //! Copy assignment operator
    ClusterRefBase& operator=(const ClusterRefBase &other) = delete;

    //! Move assignment operator
    ClusterRefBase& operator=(ClusterRefBase &&other) = default;

    const std::array<int, Level> & get_indices() const {return this->indices;}

    const int & back() const{return this->indices.back();}

  protected:
    std::array<int, Level> indices;
  private:
  };



} // rascal

#endif /* CLUSTERREFBASE_H */
