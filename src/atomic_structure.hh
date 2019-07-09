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
 * Copyright  2018  Felix Musil, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_ATOMIC_STRUCTURE_HH_
#define SRC_ATOMIC_STRUCTURE_HH_

#include "basic_types.hh"
#include "json_io.hh"

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>

// TODO(markus): CHECK for skewedness
namespace rascal {

  template <int Dim>
  struct AtomicStructure {
    /**
     * A common structure to access atom and cell related data, based on the
     * idea of the atoms object in the Atomic Simulation Environment. The
     * object contains atomic positions, the cell vectors, periodicity
     * information as well as the atomic types (element).
     *
     *  @param cell is a contiguous vector which holds the cell unit vectors.
     *
     *  @param type a vector of integers which holds the atomic type (atomic
     *  number as per periodic table).
     *
     *  @param pbc is a 0/1 vector which defines the periodicity of the given
     *  structure for each dimension
     *
     *  @param position is a vector which holds the atomic positions.
     */
    using Cell_t = Eigen::Matrix<double, Dim, Dim>;
    using Cell_ref = Eigen::Map<Cell_t>;

    using AtomTypes_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;
    using AtomTypes_ref = Eigen::Map<AtomTypes_t>;
    using ConstAtomTypes_ref = Eigen::Map<const AtomTypes_t>;

    using PBC_t = Eigen::Matrix<int, Dim, 1>;
    using PBC_ref = Eigen::Map<PBC_t>;

    using Positions_t = Eigen::Matrix<double, Dim, Eigen::Dynamic>;
    using Positions_ref = Eigen::Map<Positions_t>;

    using PositionsInput_t =
        Eigen::Ref<const Eigen::MatrixXd, 0,
                   Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

    using AtomTypesInput_t =
        Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 1>>;

    using PBCInput_t = Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 1>>;

    using ArrayB_t = Eigen::Array<bool, Eigen::Dynamic, 1>;
    using ArrayB_ref = Eigen::Ref<const ArrayB_t>;

    // Eigen types for saving atomic structure data
    Positions_t positions{};
    AtomTypes_t atom_types{};
    Cell_t cell{};
    PBC_t pbc{};
    //! Contains the information wheter an atom should be centered on or not
    ArrayB_t is_a_center_atom{};

    inline size_t get_number_of_atoms() const {
      return positions.cols();
    }

    inline void set_structure(const PositionsInput_t & positions,
                              const AtomTypesInput_t & atom_types,
                              const Eigen::Ref<const Eigen::MatrixXd> cell,
                              const PBCInput_t & pbc) {
      auto is_a_center_atom = ArrayB_t::Ones(atom_types.size());
      this->set_structure(positions, atom_types, cell, pbc, is_a_center_atom);
    }

    //! method for initializing structure data from raw Eigen types, beware:
    //! copy!
    inline void set_structure(const PositionsInput_t & positions,
                              const AtomTypesInput_t & atom_types,
                              const Eigen::Ref<const Eigen::MatrixXd> cell,
                              const PBCInput_t & pbc,
                              const ArrayB_t& is_a_center_atom) {
      // check data consistency
      auto npos{positions.cols()};
      auto ntypes{atom_types.rows()};
      auto n_center_flags{is_a_center_atom.size()};
      if (npos != ntypes or ntypes != n_center_flags) {
        std::stringstream err_str{};
        err_str << "Number of atom positions and atom types is not the same: '"
                << npos << "' != '" << ntypes << "' != '" << n_center_flags << "'.";
        throw std::runtime_error(err_str.str());
      }

      this->cell = cell;
      this->atom_types = atom_types;
      this->pbc = pbc;
      this->positions = positions;
      this->is_a_center_atom = is_a_center_atom;
    }

    // TODO(markus): add function to read from XYZ files
    inline void set_structure(const std::string & filename) {
      json j;
      std::ifstream reader(filename);
      if (not reader.is_open()) {
        throw std::runtime_error(std::string("Could not open the file: ") +
                                 filename);
      }
      reader >> j;
      reader.close();
      this->set_structure(j.begin().value());
    }

    inline void set_structure(const json & s) {
      /*
       * ASE json format is nested - here, first entry is actual data
       * structure. If in any case you should have multiple
       * <code>atoms_objects</code> in your file, which you want to read, the
       * following line has to be adapted. Nesting on the first level is
       * structure1, structure 2, etc. These could be a time series of a
       * simulation, but also just different structures you want to read in from
       * different simulation runs.  Each structure should contain the necessary
       * fields for the <code>AtomicStructure</code> object defined in the
       * header belonging to this file. Here, just the first one is read.
       */
      if (s.count("filename") == 1) {
        auto filename{s["filename"].get<std::string>()};
        this->set_structure(filename);
      } else if (s.count("cell") == 1 and s.count("atom_types") == 1 and
                 s.count("pbc") == 1 and s.count("positions") == 1) {
        auto json_atoms_object = s.get<AtomicStructure<Dim>>();
        this->set_structure(json_atoms_object);
      } else {
        throw std::runtime_error("The json input was not understood.");
      }
    }

    inline void set_structure(const AtomicStructure<Dim> & other) {
      this->positions = other.positions;
      this->atom_types = other.atom_types;
      this->cell = other.cell;
      this->pbc = other.pbc;
      this->is_a_center_atom = other.is_a_center_atom;
    }

    inline void set_structure() {}

    /**
     * Compare if another structure is identical to itself.
     *
     * Assumes that if the structure is given as json or filename related then
     * it is different. Do the comparison only if it is given as an
     * AtomicStructure or positions, pbc...
     * Used for the verlet list
     */
    inline bool is_identical(const double &) const { return true; }

    inline bool is_identical(const json_io::AtomicJsonData &, const double &) {
      return false;
    }

    inline bool is_identical(const json &, const double &) const {
      return false;
    }

    inline bool is_identical(const std::string &, const double &) const {
      return false;
    }

    inline bool is_identical(const AtomicStructure<Dim> & other,
                             const double & skin2) const {
      bool is_similar{true};
      if (this->positions.cols() == other.positions.cols()) {
        if ((this->pbc.array() != other.pbc.array()).any() or
            (this->cell.array() != other.cell.array()).any() or
            (this->is_a_center_atom != other.is_a_center_atom).any() or
            (this->positions - other.positions)
                    .rowwise()
                    .squaredNorm()
                    .maxCoeff() > skin2) {
          is_similar = false;
        }
      } else {
        is_similar = false;
      }
      return is_similar;
    }

    inline bool is_identical(const PositionsInput_t & positions,
                             const AtomTypesInput_t & /*atom_types*/,
                             const Eigen::Ref<const Eigen::MatrixXd> cell,
                             const PBCInput_t & pbc,
                             const double & skin2) const {
      bool is_similar{true};
      if (this->positions.cols() == positions.cols()) {
        if ((this->pbc.array() != pbc.array()).any() or
            (this->cell.array() != cell.array()).any() or
            /*(this->is_a_center_atom != other.is_a_center_atom).any() or*/
            (this->positions - positions).rowwise().squaredNorm().maxCoeff() >
                skin2) {
          is_similar = false;
        }
      } else {
        is_similar = false;
      }
      return is_similar;
    }
  };


  /* ---------------------------------------------------------------------- */
  template<int Dim>
  void to_json(json & j, const AtomicStructure<Dim> & s) {
    // auto cell = s.cell;
    // cell.transposeInPlace();
    j = json{{"cell", s.cell},
              {"atom_types", s.atom_types},
              {"pbc", s.pbc},
              {"positions", s.positions},
              {"is_a_center_atom", s.is_a_center_atom}};
  }

  /* ---------------------------------------------------------------------- */
  template<int Dim>
  void from_json(const json & j, AtomicStructure<Dim> & s) {
    // using Cell_t = typename AtomicStructure<Dim>::Cell_t;
    // using AtomTypes_t = typename AtomicStructure<Dim>::AtomTypes_t;
    // using PBC_t = typename AtomicStructure<Dim>::PBC_t;
    // using Positions_t = typename AtomicStructure<Dim>::Positions_t;
    using Cell_t = Eigen::MatrixXd;
    using AtomTypes_t = Eigen::VectorXi;
    using PBC_t = Eigen::VectorXi;
    using Positions_t = Eigen::MatrixXd;
    using ArrayB_t = typename AtomicStructure<Dim>::ArrayB_t;

    auto cell = j.at("cell").get<Cell_t>();
    auto atom_types = j.at("atom_types").get<AtomTypes_t>();
    auto pbc = j.at("pbc").get<PBC_t>();
    auto positions = j.at("positions").get<Positions_t>();

    cell.transposeInPlace();

    if (atom_types.size() == positions.rows()) {
      positions.transposeInPlace();
    }
    if (j.count("is_a_center_atom") == 1) {
      auto is_a_center_atom = j.at("is_a_center_atom").get<ArrayB_t>();
      s.set_structure(positions, atom_types, cell, pbc, is_a_center_atom);
    } else {
      s.set_structure(positions, atom_types, cell, pbc);
    }

  }

}  // namespace rascal

#endif  // SRC_ATOMIC_STRUCTURE_HH_
