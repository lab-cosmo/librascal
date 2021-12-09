/**
 * @file   rascal/structure_managers/cluster_ref_base.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   06 Aug 2018
 *
 * @brief implementation of non-templated base class for ClusterRef
 *
 * Copyright  2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_STRUCTURE_MANAGERS_CLUSTER_REF_BASE_HH_
#define SRC_RASCAL_STRUCTURE_MANAGERS_CLUSTER_REF_BASE_HH_

#include <cstddef>

namespace rascal {

  class ClusterRefBase {
   public:
    //! Default constructor
    ClusterRefBase() = delete;

    //! Copy constructor
    ClusterRefBase(const ClusterRefBase & other) = delete;

    //! Move constructor
    ClusterRefBase(ClusterRefBase && other) = default;

    //! Destructor
    virtual ~ClusterRefBase() = default;

    //! Copy assignment operator
    ClusterRefBase & operator=(const ClusterRefBase & other) = delete;

    //! Move assignment operator
    ClusterRefBase & operator=(ClusterRefBase && other) = delete;

    //! returns the order of the cluster
    size_t get_order() const { return this->order; }

    //! returns the layer of the cluster
    size_t get_cluster_layer() const { return this->layer; }

   protected:
    const size_t order;  //!< cluster order: atom, pair, ...
    const size_t layer;  //!< cluster layer

    //! constructor
    ClusterRefBase(size_t order, size_t layer) : order{order}, layer{layer} {}
  };

}  // namespace rascal

#endif  // SRC_RASCAL_STRUCTURE_MANAGERS_CLUSTER_REF_BASE_HH_
