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

#include <sstream>
#include <string>

namespace rascal {

  enum class SymmetryFunType { One, Gaussian, Cosine, Angular1, Angular2 };

  /* ---------------------------------------------------------------------- */
  std::string get_name(SymmetryFunType fun_type) {
    switch (fun_type) {
    case SymmetryFunType::One: {
      return "One";
      break;
    }
    case SymmetryFunType::Gaussian: {
      return "Gaussian";
      break;
    }
    case SymmetryFunType::Cosine: {
      return "Cosine";
      break;
    }
    case SymmetryFunType::Angular1: {
      return "Angular1";
      break;
    }
    case SymmetryFunType::Angular2: {
      return "Angular2";
      break;
    }
    default:
      throw std::runtime_error("undefined symmetry function type");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
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

  template<class StructureManager>
  class SymmetryFunEvaluatorBase
  {
   public:

    template <size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    //! Default constructor
    SymmetryFunEvaluatorBase() = delete;

    //! Constructor with symmetry function type
    SymmetryFunEvaluatorBase(const SymmetryFunType sym_fun_type)
        : sym_fun_type{sym_fun_type} {}

    //! Copy constructor
    SymmetryFunEvaluatorBase(const SymmetryFunEvaluatorBase &other);

    //! Move constructor
    SymmetryFunEvaluatorBase(SymmetryFunEvaluatorBase &&other) noexcept;

    //! Destructor
    virtual ~SymmetryFunEvaluatorBase() noexcept;

    //! Copy assignment operator
    SymmetryFunEvaluatorBase& operator=(const SymmetryFunEvaluatorBase &other);

    //! Move assignment operator
    SymmetryFunEvaluatorBase &
    operator=(SymmetryFunEvaluatorBase && other) noexcept;

    //! needs to be called after reading the input file and prior to the first
    //! evaluation
    virtual void init() = 0;

    //! insert a parameter (sub-)json
    void add_params(const json & params) {
      const auto & type{params.at("type").get<std::string>()};
      if (type != get_name(FunType)) {
        std::stringstream error{};
        error << "Parameter set for function type '" << type
              << "' assigned to function of type '"
              << get_name(this->sym_fun_type) << "'.";
        throw std::runtime_error(error.str());
      }
      this->raw_params.push_back(params);
    };

    //! Main worker (raison d'être)
    template <size_t Order>
    void apply(ClusterRef_t<Order> & cluster) {
      static_assert((Order == 1) or (Order == 2),
                    "only handles pairs and triplets");
      switch (Order) {
      case 1: {
        this->apply_pair(cluster);
        break;
      }
      case 2: {
        this->apply_triplet(cluster);
        break;
      }
      default: {};
      }
    }

   protected:
    virtual void apply_pair(ClusterRef_t<1> & cluster) = 0;
    virtual void apply_triplet(ClusterRef_t<2> & cluster) = 0;
    std::vector<json> raw_params{};
    const SymmetryFunType sym_fun_type;
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType FunType, class StructureManager>
  class SymmetryFunEvaluator final: public SymmetryFunEvaluatorBase {
   public:

    using Parent = SymmetryFunEvaluatorBase;
    using ParamStorage =
        Eigen::Matrix<double, FunType::NbParams, Eigen::Dynamic>;

    //! Default constructor
    SymmetryFunEvaluator() : Parent(FunType) {}

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

    void init() final;

   protected:
    static constexpr size_t AtomOrder{1};
    static constexpr size_t PairOrder{2};
    static constexpr size_t AtomLayer{
        std::get<0>(StructureManager::traits::LayerByOrder::type)};
    static constexpr size_t PairLayer{
        std::get<1>(StructureManager::traits::LayerByOrder::type)};

    ParamStorage params{};
  };


  /* ---------------------------------------------------------------------- */
  template<SymmetryFunType FunType, class StructureManager>
  void SymmetryFunEvaluator<FunType, StructureManager>::init() {
    // start by resizing the parameter storage
    this->params = params.Zero(FunType::NbParams, this->raw_params.size());

    for (const auto & raw_param: this->raw_params) {
      }
  }

}  // namespace rascal
