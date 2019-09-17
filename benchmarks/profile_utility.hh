/**
 * file  profile_scalar_cubic_spline.hh
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief Contains utilities required for the profile files.
 *
 * Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <sys/stat.h>
#include <unistd.h>
#include <fstream>

namespace rascal {
  namespace internal {

    /**
     * To write and read eigen matrix objects into binary objects. It is adapted from https://stackoverflow.com/a/25389481/10329403
     */
    template <class Matrix>
    void write_binary(const char * filename, const Matrix & matrix) {
      std::ofstream out(filename,
                        std::ios::out | std::ios::binary | std::ios::trunc);
      typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
      out.write((char *)(&rows), sizeof(typename Matrix::Index));
      out.write((char *)(&cols), sizeof(typename Matrix::Index));
      out.write((char *)matrix.data(),
                rows * cols * sizeof(typename Matrix::Scalar));
      out.close();
    }
    template <class Matrix>
    void read_binary(const char * filename, Matrix & matrix) {
      std::ifstream in(filename, std::ios::in | std::ios::binary);
      typename Matrix::Index rows = 0, cols = 0;
      in.read((char *)(&rows), sizeof(typename Matrix::Index));
      in.read((char *)(&cols), sizeof(typename Matrix::Index));
      matrix.resize(rows, cols);
      in.read((char *)matrix.data(),
              rows * cols * sizeof(typename Matrix::Scalar));
      in.close();
    }

    inline bool file_exists(const char * name) {
      struct stat buffer;
      return (stat(name, &buffer) == 0);
    }
  }  // namespace internal
}  // namespace rascal
