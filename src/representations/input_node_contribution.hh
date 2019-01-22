/**
 * file   input_node_contribution.hh
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


#ifndef _SRC_REPRESENTATIONS_INPUT_NODE_CONTRIBUTION_HH_
#define _SRC_REPRESENTATIONS_INPUT_NODE_CONTRIBUTION_HH_

#include "representations/cutoff_functions.hh"
#include "representations/symmetry_functions.hh"


namespace rascal {

  template <class StructureManager>
  class InputNodeContributionBase {
   public:

    template <size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    using StdSpecies = TupleStandardisation<StructureManager::traits::MaxOrder>;

    //! Default constructor
    InputNodeContributionBase() = delete;

    //! Constructor with symmetry function type
    InputNodeContributionBase(const SymmetryFunType sym_fun_type,
                             const size_t order)
        : sym_fun_type{sym_fun_type}, order{order} {}

    //! Copy constructor
    InputNodeContributionBase(const InputNodeContributionBase &other);

    //! Move constructor
    InputNodeContributionBase(InputNodeContributionBase &&other) noexcept;

    //! Destructor
    virtual ~InputNodeContributionBase() noexcept;

    //! Copy assignment operator
    InputNodeContributionBase& operator=(const InputNodeContributionBase &other);

    //! Move assignment operator
    InputNodeContributionBase &
    operator=(InputNodeContributionBase && other) noexcept;

    //! needs to be called after reading the input file and prior to the first
    //! evaluation. Attaches all necessary properties for precalculated values
    //! to the manager
    virtual void init(StructureManager & manager) = 0;

    //! needs to be called in the beginning of every evaluation step, refreshes
    //! the precalculated properties if necessary
    virtual void prepare(StructureManager & manager) = 0;

    //! Main worker (raison d'être)
    virtual void apply(StructureManager & manager) const = 0;

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

   protected:
    std::vector<json> raw_params{};
    const SymmetryFunType sym_fun_type;
    const size_t order;
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFun, CutoffFuntype CutFun,
            class StructureManager>
  class InputNodeContribution final
      : public InputNodeContributionBase<StructureManager> {
   public:

    using Parent = InputNodeContributionBase;
    using SymmetryFunction = SymmetryFun<SymFun>;
    using CutoffFunction = CutoffFun<CutFun>;

    //stores parameter packs ordered by cutoff radius
    using ParamStorage =
        std::map<double,
                 Eigen::Matrix<double, SymFun::NbParams, Eigen::Dynamic>>;

    constexpr static size_t MaxOrder{Parent::traits::MaxOrder};

    //! Default constructor
    InputNodeContribution() : Parent(SymFun) {}

    //! Copy constructor
    InputNodeContribution(const InputNodeContribution & other) = delete;

    //! Move constructor
    InputNodeContribution(InputNodeContribution && other) = default;

    //! Destructor
    ~InputNodeContribution() = default;

    //! Copy assignment operator
    InputNodeContribution &
    operator=(const InputNodeContribution & other) = delete;

    //! Move assignment operator
    InputNodeContribution & operator=(InputNodeContribution && other) = default;

    void init(StructureManager & manager) final;
    void prepare(StructureManager & manager) final;

    void apply(StructureManager & manager) const;

    StdSpecies get_species_combo() const; //to implement
    size_t get_index() const; //to implement

   protected:
    inline void eval_cluster(StructureManager & manager, ClusterRef_t cluster
    static constexpr size_t AtomOrder{1};
    static constexpr size_t PairOrder{2};
    static constexpr size_t AtomLayer{
        std::get<0>(StructureManager::traits::LayerByOrder::type)};
    static constexpr size_t PairLayer{
        std::get<1>(StructureManager::traits::LayerByOrder::type)};

    ParamStorage params{};
  };

}  // rascal
#include "input_node_contribution_impl.hh"

#endif  // _SRC_REPRESENTATIONS_INPUT_NODE_CONTRIBUTION_HH_
