/**
 * @file   atomic_structure.hh
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
#include "math/math_utils.hh"

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>
#include <iostream>

// TODO(markus): CHECK for skewedness
namespace rascal {

  /**
   * A common structure to access atom and cell related data, based on the
   * idea of the atoms object in the Atomic Simulation Environment. The
   * object contains atomic positions, the cell vectors, periodicity
   * information as well as the atomic types (element).
   */
  template <int Dim>
  struct AtomicStructure {
    using Cell_t = Eigen::Matrix<double, Dim, Dim>;
    using Cell_ref = Eigen::Ref<Cell_t>;
    using ConstCell_ref = const Eigen::Ref<const Cell_t>;

    using AtomTypes_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;
    using AtomTypes_ref = Eigen::Ref<AtomTypes_t>;
    using ConstAtomTypes_ref = Eigen::Ref<const AtomTypes_t>;

    using PBC_t = Eigen::Matrix<int, Dim, 1>;
    using PBC_ref = Eigen::Ref<PBC_t>;

    using Positions_t = Eigen::Matrix<double, Dim, Eigen::Dynamic>;
    using Positions_ref = Eigen::Ref<Positions_t>;

    using PositionsInput_t =
        Eigen::Ref<const Eigen::MatrixXd, 0,
                   Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

    using AtomTypesInput_t =
        Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 1>>;

    using PBCInput_t = Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 1>>;
    using CellInput_t = Eigen::Ref<const Eigen::MatrixXd>;

    using ArrayB_t = Eigen::Array<bool, Eigen::Dynamic, 1>;
    using ArrayB_ref = Eigen::Ref<const ArrayB_t>;

    template <typename T>
    using ArrayConstRef_t =
        const Eigen::Ref<const Eigen::Array<T, Eigen::Dynamic, 1>>;

    /**
     * A 3xN matrix which holds the atomic positions.
     */
    Positions_t positions{};

    /**
     * A vector of integers which holds the atomic type (atomic number as
     * per periodic table).
     */
    AtomTypes_t atom_types{};

    /**
     * A contiguous vector which holds the cell unit vectors.
     */
    Cell_t cell{};

    /**
     * A 0/1 vector which defines the periodicity of the given structure
     * for each dimension
     */
    PBC_t pbc{};

    //! Contains the information wheter an atom should be centered on or not
    //! in the form of an array of N booleans (true->center)
    ArrayB_t center_atoms_mask{};

    //! Default constructor
    AtomicStructure() = default;

    inline size_t get_number_of_atoms() const { return positions.cols(); }

    Positions_t get_scaled_positions() {
      return this->cell.inverse() * this->positions;
    }

    /**
     * fold the atoms inside the box if it has periodic boundary conditions
     */
    void wrap() {
      auto scaled_positions = this->get_scaled_positions();
      auto functor{math::MakePositivePyMod(1.)};

      for (int i_dim{0}; i_dim < Dim; ++i_dim) {
        if (this->pbc[i_dim]) {
          scaled_positions.row(i_dim) =
              scaled_positions.row(i_dim).unaryExpr(functor);
        }
      }
      this->positions = this->cell * scaled_positions;
    }

    /**
     * Set the atomic structure. The expected input are similar to the member
     * variable of the AtomicStructure class.
     *
     * By default all atoms are considered as atoms to center the
     * representation on.
     *
     * A valid atomic structure file is in the ASE json format.
     */
    inline void set_structure(const PositionsInput_t & positions,
                              const AtomTypesInput_t & atom_types,
                              const CellInput_t cell, const PBCInput_t & pbc) {
      auto center_atoms_mask = ArrayB_t::Ones(atom_types.size());
      this->set_structure(positions, atom_types, cell, pbc, center_atoms_mask);
    }

    /**
     * Method to set an property associated to the atoms like center_atoms_mask
     * to values different from the default one
     */
    template <typename T>
    inline void set_atom_property(std::string name, ArrayConstRef_t<T> array) {
      if (name == "center_atoms_mask") {
        if (atom_types.size() != array.size()) {
          throw std::runtime_error(R"(Input array does not have the same size
                                      as the number of atoms)");
        }
        center_atoms_mask = array;
      } else {
        std::string error{"The name '"};
        error += name;
        error += std::string(R"(' is not part of the possible registered
                                  fields: 'center_atoms_mask' )");
        throw std::runtime_error(error);
      }
    }

    //! method for initializing structure data from raw Eigen types, beware:
    //! copy!
    inline void set_structure(const PositionsInput_t & positions,
                              const AtomTypesInput_t & atom_types,
                              const CellInput_t cell, const PBCInput_t & pbc,
                              const ArrayB_t & center_atoms_mask) {
      // check data consistency
      auto npos{positions.cols()};
      auto ntypes{atom_types.rows()};
      auto n_center_flags{center_atoms_mask.size()};
      if (npos != ntypes or ntypes != n_center_flags) {
        std::stringstream err_str{};
        err_str << "Number of atom positions and atom types is not the same: '"
                << npos << "' != '" << ntypes << "' != '" << n_center_flags
                << "'.";
        throw std::runtime_error(err_str.str());
      }

      this->cell = cell;
      this->atom_types = atom_types;
      this->pbc = pbc;
      this->positions = positions;
      this->center_atoms_mask = center_atoms_mask;
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

      if (not s.is_object()) {
        throw std::runtime_error("The json input should be a dictionary.");
      }

      if (s.count("filename") == 1) {
        auto filename{s["filename"].get<std::string>()};
        this->set_structure(filename);
      } else if (s.count("cell") == 1 and
                 (s.count("atom_types") == 1 or s.count("numbers") == 1) and
                 s.count("pbc") == 1 and s.count("positions") == 1) {
        auto json_atoms_object = s.get<AtomicStructure<Dim>>();
        this->set_structure(json_atoms_object);
      } else {
        std::string error{
            "The json input was not understood. The input keys are: "};
        for (auto & el : s.items()) {
          error += el.key() + std::string(", ");
        }
        throw std::runtime_error(error);
      }
    }

    inline void set_structure(const AtomicStructure<Dim> & other) {
      this->positions = other.positions;
      this->atom_types = other.atom_types;
      this->cell = other.cell;
      this->pbc = other.pbc;
      this->center_atoms_mask = other.center_atoms_mask;
    }

    inline void set_structure() {}

    /**
     * Compare if another atomic structure is identical to itself.
     *
     * Assumes that if the structure is given as json or filename related then
     * it is different. Do the comparison only if it is given as an
     * AtomicStructure or positions, pbc...
     * Used for the verlet list
     *
     * @param threshold2 tolerance parameter squared for the similarity
     *                    comparison
     */
    inline bool is_similar(double) const { return true; }

    inline bool is_similar(const json_io::AtomicJsonData &, double) {
      return false;
    }

    inline bool is_similar(const json &, double) const { return false; }

    inline bool is_similar(const std::string &, double) const { return false; }

    inline bool is_similar(const AtomicStructure<Dim> & other,
                           double threshold2) const {
      bool is_similar_{true};
      if (this->positions.cols() == other.positions.cols()) {
        if ((this->pbc.array() != other.pbc.array()).any() or
            (this->cell.array() != other.cell.array()).any() or
            (this->center_atoms_mask != other.center_atoms_mask).any() or
            (this->atom_types.array() != other.atom_types.array()).any() or
            (this->positions - other.positions)
                    .rowwise()
                    .squaredNorm()
                    .maxCoeff() > threshold2) {
          is_similar_ = false;
        }
      } else {
        is_similar_ = false;
      }
      return is_similar_;
    }

    inline bool is_similar(const PositionsInput_t & positions,
                           const AtomTypesInput_t & atom_types,
                           const CellInput_t cell, const PBCInput_t & pbc,
                           double threshold2) const {
      auto center_atoms_mask = ArrayB_t::Ones(atom_types.size());
      return this->is_similar(positions, atom_types, cell, pbc,
                              center_atoms_mask, threshold2);
    }

    inline bool is_similar(const PositionsInput_t & positions,
                           const AtomTypesInput_t & atom_types,
                           const CellInput_t cell, const PBCInput_t & pbc,
                           const ArrayB_ref & center_atoms_mask,
                           double threshold2) const {
      bool is_similar_{true};
      if (this->positions.cols() == positions.cols()) {
        if ((this->pbc.array() != pbc.array()).any() or
            (this->cell.array() != cell.array()).any() or
            (this->center_atoms_mask != center_atoms_mask).any() or
            (this->atom_types.array() != atom_types.array()).any() or
            (this->positions - positions).rowwise().squaredNorm().maxCoeff() >
                threshold2) {
          is_similar_ = false;
        }
      } else {
        is_similar_ = false;
      }
      return is_similar_;
    }
  };

  /* ---------------------------------------------------------------------- */
  template <int Dim>
  void to_json(json & j, const AtomicStructure<Dim> & s) {
    auto cell = s.cell;
    cell.transposeInPlace();

    j = json{{"cell", cell},
             {"atom_types", s.atom_types},
             {"pbc", s.pbc},
             {"positions", s.positions},
             {"center_atoms_mask", s.center_atoms_mask}};
  }

  /* ---------------------------------------------------------------------- */
  template <int Dim>
  void from_json(const json & j, AtomicStructure<Dim> & s) {
    using Cell_t = Eigen::MatrixXd;
    using AtomTypes_t = Eigen::VectorXi;
    using PBC_t = Eigen::VectorXi;
    using Positions_t = Eigen::MatrixXd;
    using ArrayB_t = typename AtomicStructure<Dim>::ArrayB_t;

    auto cell = j.at("cell").get<Cell_t>();
    auto positions = j.at("positions").get<Positions_t>();

    AtomTypes_t atom_types(positions.rows());
    if (j.count("atom_types") == 1) {
      atom_types = j.at("atom_types").get<AtomTypes_t>();
    } else if (j.count("numbers") == 1) {
      atom_types = j.at("numbers").get<AtomTypes_t>();
    } else {
      throw std::runtime_error(
          R"(AtomicStructure needs atom_types or numbers keyword)");
    }
    auto pbc = j.at("pbc").get<PBC_t>();

    cell.transposeInPlace();

    if (atom_types.size() == positions.rows()) {
      positions.transposeInPlace();
    }
    if (j.count("center_atoms_mask") == 1) {
      auto center_atoms_mask = j.at("center_atoms_mask").get<ArrayB_t>();
      s.set_structure(positions, atom_types, cell, pbc, center_atoms_mask);
    } else {
      s.set_structure(positions, atom_types, cell, pbc);
    }
  }

}  // namespace rascal

#endif  // SRC_ATOMIC_STRUCTURE_HH_
