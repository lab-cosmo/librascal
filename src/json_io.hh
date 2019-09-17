/**
 * file   json_io.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   18 Jun 2018
 *
 * @brief JSON interface from nlohmanns header class
 *
 * Copyright  2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_JSON_IO_HH_
#define SRC_JSON_IO_HH_

/*
 * interface to external header-library/header-class, which makes it easy to use
 * the JSON as a first class data type. See https://github.com/nlohmann/json for
 * documentation.
 */
#include "json.hpp"
#include "rascal_utility.hh"

#include <Eigen/Dense>
#include <fstream>

// For convenience
using json = nlohmann::json;

/**
 * Utilities to convert json to Eigen types and the oposite
 */
namespace Eigen {
  namespace internal {

    template <int _Rows, int _Cols>
    struct ResizeEigen {
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_mat(
          const size_t & /* n_rows*/, const size_t & /*n_cols */,
          Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & /*s*/) {
        // assert(n_rows == _Rows);
        // assert(n_cols == _Cols);
      }
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_arr(
          const size_t & /* n_rows*/, const size_t & /*n_cols */,
          Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & /*s*/) {
        // assert(n_rows == _Rows);
        // assert(n_cols == _Cols);
      }
    };

    template <>
    struct ResizeEigen<Dynamic, Dynamic> {
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_mat(
          const size_t & n_rows, const size_t & n_cols,
          Matrix<_Scalar, Dynamic, Dynamic, _Options, _MaxRows, _MaxCols> & s) {
        s.resize(n_rows, n_cols);
      }

      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_arr(
          const size_t & n_rows, const size_t & n_cols,
          Array<_Scalar, Dynamic, Dynamic, _Options, _MaxRows, _MaxCols> & s) {
        s.resize(n_rows, n_cols);
      }
    };

    template <int _Cols>
    struct ResizeEigen<Dynamic, _Cols> {
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_mat(
          const size_t & n_rows, const size_t & /* n_cols */,
          Matrix<_Scalar, Dynamic, _Cols, _Options, _MaxRows, _MaxCols> & s) {
        // assert(n_cols == _Cols);
        s.resize(n_rows, _Cols);
      }

      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_arr(
          const size_t & n_rows, const size_t & /*n_cols */,
          Array<_Scalar, Dynamic, _Cols, _Options, _MaxRows, _MaxCols> & s) {
        // assert(n_cols == _Cols);
        s.resize(n_rows, _Cols);
      }
    };

    template <int _Rows>
    struct ResizeEigen<_Rows, Dynamic> {
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_mat(
          const size_t & /*n_rows */, const size_t & n_cols,
          Matrix<_Scalar, _Rows, Dynamic, _Options, _MaxRows, _MaxCols> & s) {
        // assert(n_rows == _Rows);
        s.resize(_Rows, n_cols);
      }

      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_arr(
          const size_t & /*n_rows */, const size_t & n_cols,
          Array<_Scalar, _Rows, Dynamic, _Options, _MaxRows, _MaxCols> & s) {
        // assert(n_rows == _Rows);
        s.resize(_Rows, n_cols);
      }
    };

    template <int _Rows, int _Cols>
    struct FillEigen {
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_mat(
          const json & j,
          Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & s) {
        for (int irow{0}; irow < s.rows(); ++irow) {
          for (int icol{0}; icol < s.cols(); ++icol) {
            s(irow, icol) = j[irow][icol].get<_Scalar>();
          }
        }
      }
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void apply_arr(
          const json & j,
          Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & s) {
        for (int irow{0}; irow < s.rows(); ++irow) {
          for (int icol{0}; icol < s.cols(); ++icol) {
            s(irow, icol) = j[irow][icol].get<_Scalar>();
          }
        }
      }
    };

    template <int _Cols>
    struct FillEigen<1, _Cols> {
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void
      apply_mat(const json & j,
                Matrix<_Scalar, 1, _Cols, _Options, _MaxRows, _MaxCols> & s) {
        for (int icol{0}; icol < s.cols(); ++icol) {
          s(icol) = j[icol].get<_Scalar>();
        }
      }
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void
      apply_arr(const json & j,
                Array<_Scalar, 1, _Cols, _Options, _MaxRows, _MaxCols> & s) {
        for (int icol{0}; icol < s.cols(); ++icol) {
          s(icol) = j[icol].get<_Scalar>();
        }
      }
    };

    template <int _Rows>
    struct FillEigen<_Rows, 1> {
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void
      apply_mat(const json & j,
                Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, _MaxCols> & s) {
        for (int icol{0}; icol < s.rows(); ++icol) {
          s(icol) = j[icol].get<_Scalar>();
        }
      }
      template <typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
      static void
      apply_arr(const json & j,
                Array<_Scalar, _Rows, 1, _Options, _MaxRows, _MaxCols> & s) {
        for (int icol{0}; icol < s.rows(); ++icol) {
          s(icol) = j[icol].get<_Scalar>();
        }
      }
    };

  }  // namespace internal
  /* ---------------------------------------------------------------------- */
  /**
   * By convention 1D array are stored as a single vector and not a matrix
   * with one dimension set to 1.
   */
  template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows,
            int _MaxCols>
  void to_json(
      json & j,
      const Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & s) {
    if (_Rows != 1 and _Cols != 1) {
      for (int irow{0}; irow < s.rows(); ++irow) {
        j.push_back(json::array());
        for (int icol{0}; icol < s.cols(); ++icol) {
          j[irow][icol] = s(irow, icol);
        }
      }
    } else if (_Rows == 1) {
      for (int icol{0}; icol < s.cols(); ++icol) {
        j[icol] = s(icol);
      }
    } else if (_Cols == 1) {
      for (int irow{0}; irow < s.rows(); ++irow) {
        j[irow] = s(irow);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows,
            int _MaxCols>
  void
  from_json(const json & j,
            Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & s) {
    if (j.is_array() != true) {
      std::string error{"the input json is not an array "};
      std::string error2{j.dump(2)};
      throw std::runtime_error(error + error2);
    }
    size_t n_rows{j.size()};
    size_t n_cols{j[0].size()};
    for (size_t irow{0}; irow < n_rows; ++irow) {
      if (n_cols != j[irow].size()) {
        std::string error{"the input json does not have cols of same sizes"};
        throw std::runtime_error(error);
      }
    }

    internal::ResizeEigen<_Rows, _Cols>::apply_mat(n_rows, n_cols, s);

    internal::FillEigen<_Rows, _Cols>::apply_mat(j, s);
  }

  /* ---------------------------------------------------------------------- */
  template <typename _Scalar, int _Rows, int _Cols, int _Options>
  void to_json(json & j, const Array<_Scalar, _Rows, _Cols, _Options> & s) {
    if (_Rows != 1 and _Cols != 1) {
      for (int irow{0}; irow < s.rows(); ++irow) {
        j.push_back(json::array());
        for (int icol{0}; icol < s.cols(); ++icol) {
          j[irow][icol] = s(irow, icol);
        }
      }
    } else if (_Rows == 1) {
      for (int icol{0}; icol < s.cols(); ++icol) {
        j[icol] = s(icol);
      }
    } else if (_Cols == 1) {
      for (int irow{0}; irow < s.rows(); ++irow) {
        j[irow] = s(irow);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  template <typename _Scalar, int _Rows, int _Cols, int _Options>
  void from_json(const json & j, Array<_Scalar, _Rows, _Cols, _Options> & s) {
    if (j.is_array() != true) {
      std::string error{"the input json is not an array "};
      std::string error2{j.dump(2)};
      throw std::runtime_error(error + error2);
    }
    size_t n_rows{j.size()};
    size_t n_cols{j[0].size()};
    for (size_t irow{0}; irow < n_rows; ++irow) {
      if (n_cols != j[irow].size()) {
        std::string error{"the input json does not have cols of same sizes"};
        throw std::runtime_error(error);
      }
    }

    internal::ResizeEigen<_Rows, _Cols>::apply_arr(n_rows, n_cols, s);

    internal::FillEigen<_Rows, _Cols>::apply_arr(j, s);
  }
}  // namespace Eigen

/*
 * All functions and classes are in the namespace <code>rascal</code>, which
 * ensures that they don't clash with other libraries one might use in
 * conjunction.
 */
namespace rascal {
  namespace json_io {

    //! load a json file
    json load(const std::string & filename);

    //! load a json file in text format
    json load_txt(const std::string & filename);

    //! load a json file in ubjson binary format
    json load_bin(const std::string & filename);

    /**
     * Object to deserialize the content of a JSON file containing Atomic
     * Simulation Environment (ASE) type atomic structures, the nlohmann::json
     * needs a <code>struct</code> with standard data types for deserialization
     */
    struct AtomicJsonData {
      /**
       *  @param cell is a vector a vector of vectors which holds the cell unit
       *  vectors.
       *
       *  @param type a vector of integers which holds the atomic type (atomic
       *  number from periodic table).
       *
       *  @param pbc is a 0/1 vector which says, where periodic boundary
       *  conditions are applied.
       *
       *  @param position is a vector of vectors which holds the atomic
       *  positions.
       */
      std::vector<std::vector<double>> cell{};
      std::vector<int> type{};
      std::vector<int> pbc{};
      std::vector<std::vector<double>> position{};
    };

    /**
     * Function to convert to a JSON object format with the given keywords. It
     * is an overload of the function defined in the header class
     * json.hpp.
     */
    void to_json(json & j, AtomicJsonData & s);

    /**
     * Function used to read from the JSON file, given the keywords and convert
     * the data into standard types. Overload of the function defined in
     * json.hpp class header.
     */
    void from_json(const json & j, AtomicJsonData & s);

    /**
     * checks a value-unit pair of form {"value": 5.6, "unit": "Å"} against the
     * expected unit, returns value if successful and throws an error if not
     */
    double check_units(const std::string & expected_unit,
                       const json & parameter);
  }  // namespace json_io
}  // namespace rascal

#endif  // SRC_JSON_IO_HH_
