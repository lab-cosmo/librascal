/**
 * file   representation_manager_behler_parinello_impl.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   13 Dec 2018
 *
 * @brief  implementation for Behler-Parinello representation manager
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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


#include "utils/for_each_at_order.hh"
namespace rascal {

  template <class StructureManager>
  BehlerParinello<StructureManager>::BehlerParinello(
      StructureManager & structure, const json & hypers)
      : structure{structure}, species{structure} {}

  
  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  void BehlerParinello<StructureManager>::compute() {
    using utils::for_each_at_order;

    for (const auto && species_key_val : this->symmetry_functions) {
      const auto & species_combo{key_val.first()};
      const auto & functions_by_cutoff{key_val.second()};

      auto & manager{this->species[species_combo]};
      const auto & order{species_combo.get_order()};

      switch (order) {
      case 2: {
        for_each_at_order<2>(function, manager, functions_by_cutoff);
        break;
      }
      case 3: {
        for_each_at_order<3>(function, manager, functions_by_cutoff);
        break;
      }
      }
    }
  }

}  // namespace rascal
