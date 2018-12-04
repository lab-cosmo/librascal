/**
 * file   atomic_structure.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   08 August 2018
 *
 * @brief common data type for atomic structure data including positions, types,
 *        cell and periodic boundary conditions
 *
 * Copyright © 2018  Felix Musil, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef ATOMIC_STRUCTURE_H
#define ATOMIC_STRUCTURE_H

#include "basic_types.hh"
#include "json_io.hh"

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>
#include <cmath>




// TODO(markus): CHECK for skewedness and add a function for rescaling the
// positions for internal use
namespace rascal {

  template <int Dim>
  struct AtomicStructure {
    /**
     * A common structure to access atom and cell related data, based on the
     * idea of the atoms object in the Atomic Simulation Environment. The object
     * contains atomic positions, the cell vectors, periodicity information as
     * well as the atomic types (element).
     *
     *  \param cell is a contiguous vector which holds the cell unit vectors.
     *
     *  \param type a vector of integers which holds the atomic type (atomic
     *  number as per periodic table).
     *
     *  \param pbc is a 0/1 vector which defines the periodicity of the given
     *  structure for each dimension
     *
     *  \param position is a vector which holds the atomic positions.
     */
    using Cell_t = Eigen::Matrix<double, Dim, Dim>;
    using Cell_ref = Eigen::Map<Cell_t>;

    using AtomTypes_t = Eigen::Matrix<int, 1, Eigen::Dynamic>;
    using AtomTypes_ref = Eigen::Map<AtomTypes_t>;
    using ConstAtomTypes_ref = Eigen::Map<const AtomTypes_t>;

    using PBC_t = Eigen::Matrix<int, Dim, 1>;
    using PBC_ref = Eigen::Map<PBC_t>;

    using Positions_t = Eigen::Matrix<double, Dim,
                                      Eigen::Dynamic>;
    using Positions_ref = Eigen::Map<Positions_t>;

    using PositionsInput_t =
      Eigen::Ref<const Eigen::MatrixXd, 0,
                 Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

    using AtomTypesInput_t =
      Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 1>>;

    using PBCInput_t =
      Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 1>>;

    // Eigen types for saving atomic structure data
    Positions_t positions;
    AtomTypes_t atom_types;
    Cell_t cell;
    PBC_t pbc;

    //! method for initializing structure data from raw Eigen types, beware:
    //! copy!
    void set_structure(const PositionsInput_t & positions,
                       const AtomTypesInput_t &  atom_types,
                       const Eigen::Ref<const Eigen::MatrixXd> cell,
                       const PBCInput_t & pbc) {
      // check data consistency
      auto npos{positions.cols()};
      auto ntypes{atom_types.rows()};
      if (npos != ntypes) {
        std::stringstream err_str{};
        err_str << "Number of atom positions and atom types is not the same: '"
                << npos
                << "' != '" << ntypes << "'.";
        throw std::runtime_error(err_str.str());
      }

      this->cell = cell;
      this->atom_types = atom_types;
      this->pbc = pbc;
      this->positions = positions;
    }

    //! method for initializing structure from a json object; data is copied
    void set_structure(const json_io::AtomicJsonData & s) {
      // internal std::vector for reading from json, necessary for push_back, no
      // direct mapping possible
      std::vector<double> cell_data{};
      std::vector<int> type_data{};
      std::vector<int> pbc_data{};
      std::vector<double> pos_data{};

      // check for empty data set
      try {
        auto pos_size = s.position.size();
        if (pos_size == 0) {
          throw std::runtime_error("No atomic structure defined. "
                                   "Read structure first!");
        }
      } catch (const std::exception & e) {
        std::cerr << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // get data out of the json object to access with Eigen::Map
      for (auto vec : s.cell) {
        for (auto coord : vec) {
          cell_data.push_back(coord);
        }
      }
      // elements
      for (auto val : s.type) {
        type_data.push_back(val);
      }
      // periodicity
      for (auto val : s.pbc) {
        pbc_data.push_back(val);
      }
      // positions
      for (auto pos : s.position) {
        for (auto coord : pos) {
          pos_data.push_back(coord);
        }
      }

      // check data consistency
      auto npos{positions.size()/Dim};
      auto ntypes{atom_types.size()};
      if (npos != ntypes) {
        std::stringstream err_str{};
        err_str << "Number of atom positions and atom types is not the same: '"
                << npos
                << "' != '" << ntypes << "'.";
        throw std::runtime_error(err_str.str());
      }
      // associate them to internal data structure
      this->cell = Cell_ref(cell_data.data());
      this->atom_types = AtomTypes_ref(type_data.data(), type_data.size());
      this->pbc = PBC_ref(pbc_data.data());
      this->positions = Positions_ref(pos_data.data(), Dim,
                                      pos_data.size() / Dim);
    }
  };
} // rascal

#endif /* ATOMIC_STRUCTURE_H */
