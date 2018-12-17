/**
 * file   symmetry_functions.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   17 Dec 2018
 *
 * @brief implementation of symmetry functions for neural nets (G-functions in
 * Behler-Parinello-speak)
 *
 * Copyright © 2018 Till Junge, Markus Stricker COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/property_typed.hh"

#include "Eigen/Dense"
#include "json_io.hh"

namespace rascal {

  enum class SymmetryFunType { One, Gaussian, Cosine, Angular1, Angular2 };

  template <SymmetryFunType FunType>
  struct SymmetryFun {};

  template <>
  struct SymmetryFun<SymmetryFunType::Gaussian> {

    static constexpr size_t NbParams{2};
    using ParamShape = Eigen::MatrixBase<double, NbParams, 1>;
    /**
     * usually, derivatives are alligned with the distance vector, in which case
     * a scalar return type is sufficient. (important for triplet-related
     * functions)
     */
    static constexpr bool DerivativeIsCollinear{true};

    static double eval_function(const Eigen::MatrixBase<ParamShape> & params,
                                const double & r_ij) {
      auto && eta{params(0)};
      auto && r_s{params(1)};
      auto && delta_r = r_ij - r_s;
      return exp(-eta * delta_r * delta_r);
    }

    template <class Derived>
    static auto eval_derivative(const Eigen::MatrixBase<ParamShape> & params,
                                const double & r_ij,
                                const Eigen::MatrixBase<Derived> & n_ij)
        -> decltype(auto) {
      auto && eta{params(0)};
      auto && r_s{params(1)};
      auto && delta_r{r_ij - r_s};
      return n_ij * (-2. * eta * delta_r * exp(-eta * delta_r * delta_r));
    }
  };

  template <SymmetryFunType FunType, class StructureManager>
  class SymmetryFunEvaluator final {
   public:
    using ParamStorage =
        Eigen::Matrix<double, FunType::NbParams, Eigen::Dynamic>;

    //! Default constructor
    SymmetryFunEvaluator() = delete;

    //! construction from json input
    SymmetryFunEvaluator(const json & hypers);

    //! Copy constructor
    SymmetryFunEvaluator(const SymmetryFunEvaluator & other) = delete;

    //! Move constructor
    SymmetryFunEvaluator(SymmetryFunEvaluator && other) = default;

    //! Destructor
    ~SymmetryFunEvaluator() = default;

    //! Copy assignment operator
    SymmetryFunEvaluator &
    operator=(const SymmetryFunEvaluator & other) = delete;

    //! Move assignment operator
    SymmetryFunEvaluator & operator=(SymmetryFunEvaluator && other) = default;

   protected:
    static constexpr size_t AtomOrder{1};
    static constexpr size_t PairOrder{2};
    static constexpr size_t AtomLayer{
        std::get<0>(StructureManager::traits::LayerByOrder::type)};
    static constexpr size_t PairLayer{
        std::get<1>(StructureManager::traits::LayerByOrder::type)};

    ParamStorage params;
    TypedProperty<double, AtomOrder, AtomLayer> function_values;
    TypedProperty<double, PairOrder, PairLayer> derivative_values;

   private:
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Constructor to size the property according to json objects content
   */
  template <SymmetryFunType FunType, class StructureManager>
  SymmetryFunEvaluator<FunType, StructureManager>::SymmetryFunEvaluator(
      const json & hypers)
      : {
    // get number of hyper parameters from json object

    // Nested structure for similar hyper parameters? Is `hyper` the complete
    // json input file?

    // NbParams is already sized?
    // NbParams is different, with pair or triplets

    std::vector<std::vector<double>> params_tmp{
        hypers.at("params").get<std::vector<std::vector<double>>>()};
    size_t n_symfun{params_tmp.size()};

    // possible to directly construct `params` from json or is this intermediate
    // step necessary?
  }

}  // namespace rascal
