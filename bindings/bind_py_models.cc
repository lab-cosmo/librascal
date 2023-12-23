/**
 * @file   bind_py_models.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   16 Jun 2019
 *
 * @brief  File for binding the models
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "bind_py_models.hh"

namespace rascal {

  //! Register kernel classes
  template <class Kernel>
  static py::class_<Kernel> add_kernel(py::module & mod,
                                       py::module & /*m_internal*/) {
    std::string kernel_name = internal::GetBindingTypeName<Kernel>();
    py::class_<Kernel> kernel(mod, kernel_name.c_str());
    kernel.def(py::init([](const py::dict & hyper) {
      // convert to json
      json hypers = hyper;
      return std::make_unique<Kernel>(hypers);
    }));

    return kernel;
  }

  //! Register compute functions of the Kernel class
  template <internal::KernelType Type, class Calculator,
            class StructureManagers, class CalculatorBind>
  void bind_kernel_compute_function(CalculatorBind & kernel) {
    kernel.def("compute",
               py::overload_cast<const Calculator &, const StructureManagers &,
                                 const StructureManagers &>(
                   &Kernel::template compute<Calculator, StructureManagers>),
               py::call_guard<py::gil_scoped_release>(),
               R"(Compute the kernel between two sets of atomic structures,
              i.e. StructureManagerCollections. The representation of the
              atomic structures computed with calculator should have already
              been computed.)");
    kernel.def("compute",
               py::overload_cast<const Calculator &, const StructureManagers &>(
                   &Kernel::template compute<Calculator, StructureManagers>),
               py::call_guard<py::gil_scoped_release>(),
               R"(Compute the kernel between a set of atomic structures,
              i.e. StructureManagerCollections, and itself. The representation
              of the atomic structures computed with calculator should have
              already been computed.)");
  }

  //! Register compute functions of the SparseKernel class
  template <internal::SparseKernelType Type, class Calculator,
            class StructureManagers, class SparsePoints, class CalculatorBind>
  void bind_sparse_kernel_compute_function(CalculatorBind & kernel) {
    kernel.def(
        "compute",
        py::overload_cast<const Calculator &, const StructureManagers &,
                          const SparsePoints &>(
            &SparseKernel::template compute<Calculator, StructureManagers,
                                            SparsePoints>),
        py::call_guard<py::gil_scoped_release>(),
        R"(Compute the sparse kernel between the representation of a set of
            atomic structures, i.e. StructureManagerCollections, and a set of
            SparsePoints, i.e. the basis used by the sparse method.
            The representation of the atomic structures computed with Calculator
            should have already been computed.)");
    kernel.def(
        "compute_local",
        py::overload_cast<const Calculator &, const StructureManagers &,
                          const SparsePoints &>(
            &SparseKernel::template compute_local<Calculator, StructureManagers,
                                            SparsePoints>),
        py::call_guard<py::gil_scoped_release>(),
        R"(Compute the sparse kernel between the representation of a set of
            atomic structures, i.e. StructureManagerCollections, and a set of
            SparsePoints, i.e. the basis used by the sparse method.
            The representation of the atomic structures computed with Calculator
            should have already been computed.)");
    kernel.def("compute",
               py::overload_cast<const SparsePoints &>(
                   &SparseKernel::template compute<SparsePoints>),
               py::call_guard<py::gil_scoped_release>(),
               R"(Compute the kernel between a set of SparsePoints, i.e.
            the basis used by the sparse method, and itself.)");
    kernel.def(
        "compute_derivative",
        &SparseKernel::template compute_derivative<
            Calculator, StructureManagers, SparsePoints>,
        py::call_guard<py::gil_scoped_release>(),
        R"(Compute the sparse kernel between the gradient of representation of a
            set of atomic structures w.r.t. the atomic positions,
            i.e. StructureManagerCollections, and a set of SparsePoints, i.e.
            the basis used by the sparse method. The gradients of the
            representation of the atomic structures computed with Calculator
            should have already been computed.)");
  }

  //! Register a pseudo points class
  template <class SparsePoints>
  static py::class_<SparsePoints>
  add_sparse_points(py::module & mod, py::module & /*m_internal*/) {
    std::string pseudo_point_name =
        internal::GetBindingTypeName<SparsePoints>();
    py::class_<SparsePoints> sparse_points(mod, pseudo_point_name.c_str());
    sparse_points.def(py::init());
    sparse_points.def("size", &SparsePoints::size);
    sparse_points.def("get_features", &SparsePoints::get_features);
    return sparse_points;
  }

  //! register the population mechanism of the pseudo points class
  template <class ManagerCollection, class Calculator, class SparsePoints>
  void bind_sparse_points_push_back(py::class_<SparsePoints> & sparse_points) {
    using ManagerPtr_t = typename ManagerCollection::value_type;
    using Manager_t = typename ManagerPtr_t::element_type;
    sparse_points.def(
        "extend",
        py::overload_cast<const Calculator &, const ManagerCollection &,
                          const std::vector<std::vector<int>> &>(
            &SparsePoints::template push_back<ManagerCollection>),
        py::call_guard<py::gil_scoped_release>());
    sparse_points.def(
        "extend",
        py::overload_cast<const Calculator &, std::shared_ptr<Manager_t>,
                          const std::vector<int> &>(
            &SparsePoints::template push_back<Manager_t>),
        py::call_guard<py::gil_scoped_release>());
  }

  //! register the population mechanism of the pseudo points class
  template <class KernelImpl, class Calculator, class Managers,
            class SparsePoints>
  void bind_compute_numerical_kernel_gradients(py::module & mod) {
    mod.def("compute_numerical_kernel_gradients",
            &compute_numerical_kernel_gradients<KernelImpl, Calculator,
                                                Managers, SparsePoints>,
            py::call_guard<py::gil_scoped_release>());
  }

  template <class ManagerCollection, class Calculator, class SparsePoints>
  void bind_compute_gradients(py::module & mod, py::module & /*m_internal*/) {
    using Manager_t = typename ManagerCollection::Manager_t;
    mod.def(
        "compute_sparse_kernel_gradients",
        [](const Calculator & calculator, SparseKernel & kernel,
           ManagerCollection & managers, SparsePoints & sparse_points,
           math::Vector_t & weights) {
          std::string force_name = compute_sparse_kernel_gradients(
              calculator, kernel, managers, sparse_points, weights);
          size_t n_centers{0};
          // find the total number of gradients
          for (const auto & manager : managers) {
            n_centers += manager->size();
          }
          math::Matrix_t gradients_global{n_centers, ThreeD};
          size_t i_center{0};
          for (const auto & manager : managers) {
            auto && gradients{*manager->template get_property<
                Property<double, 1, Manager_t, 1, ThreeD>>(force_name, true)};
            gradients_global.block(i_center, 0, manager->size(), ThreeD) =
                gradients.view();
            i_center += manager->size();
          }
          return gradients_global;
        },
        py::call_guard<py::gil_scoped_release>());

    mod.def(
        "compute_sparse_kernel_neg_stress",
        [](const Calculator & calculator, SparseKernel & kernel,
           ManagerCollection & managers, SparsePoints & sparse_points,
           math::Vector_t & weights) {
          std::string neg_stress_name = compute_sparse_kernel_neg_stress(
              calculator, kernel, managers, sparse_points, weights);
          math::Matrix_t neg_stress_global{managers.size(), 6};
          size_t i_manager{0};
          for (const auto & manager : managers) {
            auto && neg_stress{
                *manager
                     ->template get_property<Property<double, 0, Manager_t, 6>>(
                         neg_stress_name, true)};
            neg_stress_global.block(i_manager, 0, 1, 6) =
                Eigen::Map<const math::Matrix_t>(neg_stress.view().data(), 1,
                                                 6);
            i_manager++;
          }
          return neg_stress_global;
        },
        py::call_guard<py::gil_scoped_release>());

    mod.def(
        "compute_sparse_kernel_local_neg_stress",
        [](const Calculator & calculator, SparseKernel & kernel,
           ManagerCollection & managers, SparsePoints & sparse_points,
           math::Vector_t & weights) {
          std::string neg_stress_name = compute_sparse_kernel_local_neg_stress(
              calculator, kernel, managers, sparse_points, weights);
          size_t nb_centers_over_managers{0};
          for (const auto & manager : managers) {
            
            nb_centers_over_managers+=manager->size();
            //for (auto center : manager) {
            //  nb_centers_over_managers++;
            //}
          }
          math::Matrix_t local_neg_stress{nb_centers_over_managers, 9};

          size_t i_center{0};
          size_t i_manager_nb_center{0};
          for (const auto & manager : managers) {
            i_manager_nb_center=manager->size();
            // i_manager_nb_center = 0;
            // for (auto center : manager) {
            //   i_manager_nb_center++;
            // }
            auto && neg_stress{
                *manager
                    ->template get_property<Property<double, 1, Manager_t, 9>>(
                        neg_stress_name, true, true, true)};
            local_neg_stress.block(i_center, 0, i_manager_nb_center, 9) =
                Eigen::Map<const math::Matrix_t>(neg_stress.view().data(), 1, 9);
            i_center += i_manager_nb_center;
          }
          return local_neg_stress;
        },
        py::call_guard<py::gil_scoped_release>());
  }

  /**
   * Function to bind the representation managers to python
   *
   * @param mod pybind11 representation of the python module the represenation
   *             managers will be included to
   * @param m_internal pybind11 representation of the python module that are
   *                   needed but not useful to use on the python side
   *
   */
  void add_models(py::module & mod, py::module & m_internal) {
    py::module m_kernels = mod.def_submodule("kernels");
    m_kernels.doc() = "Collection of Kernels";

    // Defines a particular structure manager type
    using ManagerCollection_1_t =
        ManagerCollection<StructureManagerCenters, AdaptorNeighbourList,
                          AdaptorStrict>;
    using ManagerCollection_2_t =
        ManagerCollection<StructureManagerCenters, AdaptorNeighbourList,
                          AdaptorCenterContribution, AdaptorStrict>;
    // Defines the representation manager type for the particular structure
    // manager
    using Calc1_t = CalculatorSphericalInvariants;
    using SparsePoints_1_t = SparsePointsBlockSparse<Calc1_t>;

    // Bind the interface of this representation manager
    auto kernel = add_kernel<Kernel>(m_kernels, m_internal);
    internal::bind_dict_representation(kernel);
    bind_kernel_compute_function<internal::KernelType::Cosine, Calc1_t,
                                 ManagerCollection_1_t>(kernel);
    bind_kernel_compute_function<internal::KernelType::Cosine, Calc1_t,
                                 ManagerCollection_2_t>(kernel);

    // bind the sparse kernel and pseudo points class
    auto sparse_kernel = add_kernel<SparseKernel>(m_kernels, m_internal);
    internal::bind_dict_representation(sparse_kernel);
    bind_sparse_kernel_compute_function<internal::SparseKernelType::GAP,
                                        Calc1_t, ManagerCollection_2_t,
                                        SparsePoints_1_t>(sparse_kernel);
    auto sparse_points =
        add_sparse_points<SparsePoints_1_t>(m_kernels, m_internal);
    bind_sparse_points_push_back<ManagerCollection_2_t, Calc1_t>(sparse_points);
    internal::bind_dict_representation(sparse_points);

    bind_compute_gradients<ManagerCollection_2_t, Calc1_t, SparsePoints_1_t>(
        mod, m_internal);
    bind_compute_numerical_kernel_gradients<
        SparseKernel, Calc1_t, ManagerCollection_2_t, SparsePoints_1_t>(mod);
  }
}  // namespace rascal
