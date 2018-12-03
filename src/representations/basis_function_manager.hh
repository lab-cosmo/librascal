/**
 * file   basis_function_manager.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   08 May 2018
 *
 * @brief manager for basis functions, which are the input of the neural
 *        network. e.g. symmetry functions or spherical harmonics
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef BASIS_FUNCTION_MANAGER_H
#define BASIS_FUNCTION_MANAGER_H

namespace rascal {
  class BasisFunManager {
   public:
    enum class BasisFunType: int;
    enum class CutoffFunType: int;
    using uint = insigned int;

    //! Default constructor
    BasisFunManager() = delete;

    //! Construct from file
    explicit BasisFunManager(FILE *);

    //! Copy constructor
    BasisFunManager(const BasisFunManager &other) = delete;

    //! Move constructor
    BasisFunManager(BasisFunManager &&other) = default;

    //! Destructor
    virtual ~BasisFunManager()  = default;

    //! Copy assignment operator
    BasisFunManager& operator=(const BasisFunManager &other) = delete;

    //! Move assignment operator
    BasisFunManager& operator=(BasisFunManager && other) = default;

    //! Initialize hyperparameters
    void read_hyperparameters(const FILE *);
    void set_random_hyperparameters(const int);

    //! Number of hyperparameters per basis function type
    constexpr static unit get_nhyper(const BasisFunType& fun_type);

    //! Basis functions and derivative
    template<BasisFunType fun_type>
    inline double comp_fun(const double * const param, const double * rij);
    //! return a matrix
    template<BasisFunType func_type, T>
    decltype(auto) comp_Dfun(const double * const param,
           const double * const rij);
// inline double comp_Dfun(const double * const param,
//     const double * const rij);
  }

}  // namespace rascal

#endif /* BASIS_FUNCTION_MANAGER_H */
