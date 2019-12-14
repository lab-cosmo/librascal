/**
 * file   behler_feature.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief utilities for the definition of input nodes (combinations of symmetry
 * functions with cutoff functions
 *
 * Copyright © 2019 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_HH_
#define SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_HH_

#include "rascal/json_io.hh"
#include "rascal/representations/cutoff_functions.hh"
#include "rascal/representations/symmetry_functions.hh"
#include "rascal/structure_managers/property.hh"

namespace rascal {

  template <SymmetryFunType... SymFunTypes, CutoffFunctionType... CutFunTypes>
  class BehlerFeatureBase {
   public:
    constexpr static int MaxBehlerOrder{3};
    enum class RepeatedSpecies {
      Unknown,    // has not been evaluated yet
      Not,        // all species in this cluster are unique
      All,        // all atoms in this cluster are same species
      FirstTwo,   // the first two atoms of this cluster are of same species
      SecondTwo,  // the second two atoms of this cluster are of same species
      OuterTwo
    };  // the first and last atom in this cluster are of same species
    using StdSpecies = TupleStandardisation<int, MaxBehlerOrder>;

    using Hypers_t = json;

    //! Default constructor
    BehlerFeatureBase() = delete;

    //! Constructor with symmetry function type
    BehlerFeatureBase(const SymmetryFunType & sym_fun_type,
                      const internal::CutoffFunctionType & cut_fun_type,
                      const size_t & order, const Hypers_t & raw_params)
        : sym_fun_type{sym_fun_type}, cut_fun_type{cut_fun_type}, order{order},
          raw_params{raw_params}, output_property{output_property} {}

    //! Copy constructor
    BehlerFeatureBase(const BehlerFeatureBase & other);

    //! Move constructor
    BehlerFeatureBase(BehlerFeatureBase && other) noexcept;

    //! Destructor
    virtual ~BehlerFeatureBase() noexcept;

    //! Copy assignment operator
    BehlerFeatureBase & operator=(const BehlerFeatureBase & other);

    //! Move assignment operator
    BehlerFeatureBase & operator=(BehlerFeatureBase && other) noexcept;

    //! needs to be called after reading the input file and prior to the first
    //! evaluation. Attaches all necessary properties for precalculated values
    //! to the manager
    virtual void init(const UnitStyle & units) = 0;

    //! Main worker (raison d'être) computes input node values
    template <class StructureManager>
    inline void compute(StructureManager & manager,
                        std::shared_ptr<PropertyBase> output_values) const;

    //! Main worker (raison d'être) computes input node values and derivatives
    template <class StructureManager>
    inline void compute(StructureManager & manager,
                        std::shared_ptr<PropertyBase> output_values,
                        std::shared_ptr<PropertyBase> output_derivatives) const;

    //! insert a parameter (sub-)json
    void add_params(const json & params) {
      const auto & type{params.at("type").get<std::string>()};
      if (type != get_name(this->sym_fun_type)) {
        std::stringstream error{};
        error << "Parameter set for function type '" << type
              << "' assigned to function of type '"
              << get_name(this->sym_fun_type) << "'.";
        throw std::runtime_error(error.str());
      }
      this->raw_params.push_back(params);
    }

   protected:
    template <SymmetryFunType... FunTypes_>
    class SymFunctionsVTable;

    template <SymmetryFunType SymFunType,
              internal::CutoffFunctionType... CutoffFunTypes_>
    class CutoffFunctionsVTable;

    template <SymmetryFunType SymFunType, class StructureManager>
    void compute_helper(StructureManager & manager,
                        std::shared_ptr<PropertyBase> output_values) const;

    template <SymmetryFunType SymFunType, class StructureManager>
    void compute_helper(StructureManager & manager,
                        std::shared_ptr<PropertyBase> output_values,
                        std::shared_ptr<PropertyBase> output_derivatives) const;

    const SymmetryFunType sym_fun_type;
    const internal::CutoffFunctionType cut_fun_type;
    const size_t order;
    std::vector<json> raw_params{};
    RepeatedSpecies species_repetition{RepeatedSpecies::Unknown};
    StdSpecies species_combo{};
    bool is_initialised{false};
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, internal::CutoffFunctionType CutFunType>
  class BehlerFeature final : public BehlerFeatureBase {
   public:
    using Parent = BehlerFeatureBase;
    using SymmetryFunction = SymmetryFun<SymFunType>;
    using CutoffFunction = internal::CutoffFunction<CutFunType>;
    using StdSpecies = typename Parent::StdSpecies;

    // stores parameter packs ordered by cutoff radius
    using ParamStorage =
        std::map<double,
                 Eigen::Matrix<double, SymmetryFun<SymFunType>::NbParams,
                               Eigen::Dynamic>>;

    //! Default constructor
    BehlerFeature() : Parent(SymFunType, CutFunType, SymmetryFunction::Order) {}

    //! Copy constructor
    BehlerFeature(const BehlerFeature & other) = delete;

    //! Move constructor
    BehlerFeature(BehlerFeature && other) = default;
    //! Destructor
    ~BehlerFeature() = default;

    //! Copy assignment operator
    BehlerFeature & operator=(const BehlerFeature & other) = delete;

    //! Move assignment operator
    BehlerFeature & operator=(BehlerFeature && other) = default;

    void init(const UnitStyle & units) final;

    template <class StructureManager>
    void compute(StructureManager & manager,
                 std::shared_ptr<PropertyBase> output) const;

    template <class StructureManager>
    void compute(StructureManager & manager,
                 std::shared_ptr<PropertyBase> output,
                 std::shared_ptr<PropertyBase> output_derivatives) const;

    size_t get_index() const;  // to implement

   protected:
    template <class StructureManager, size_t Order>
    inline double
    eval_cluster(StructureManager & manager,
                 const typename StructureManager::template ClusterRef_t<Order> &
                     cluster);
    ParamStorage params{};
    CutoffFunction cuf_off_fun;
  };

}  // namespace rascal
#include "behler_feature_impl.hh"

#endif  // SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_HH_
