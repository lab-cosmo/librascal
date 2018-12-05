/**
 * @file   bind_py_neighbour_manager.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   9 Mai 2018
 *
 * @brief  File for binding the Neighbour Managers
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#include "bind_include.hh"

template<typename SMI>
using LayerByOrder = typename SMI::traits::LayerByOrder;

template<typename SMI, size_t Order>
struct HelperLayer {
  static constexpr size_t layer{
              compute_cluster_layer<Order>(LayerByOrder<SMI>{})};
};

template<typename SMI, size_t Order>
using ClusterRefKey_t = ClusterRefKey<Order, HelperLayer<SMI, Order>::layer>;

template<typename SMI, size_t Order>
using ClusterRef_t =
    typename StructureManager<SMI>::template ClusterRef<Order>;


template<typename SMI, size_t Order>
using PyClusterRef = py::class_<
              ClusterRef_t<SMI, Order>, ClusterRefKey_t<SMI, Order>>;

template<typename StructureManagerImplementation>
using PyManager = py::class_<StructureManagerImplementation,
      StructureManager<StructureManagerImplementation>>;


template<size_t Order, size_t Layer>
void add_cluster_ref(py::module & m) {
  std::string cluster_parent_name =
      internal::GetBindingTypeName<ClusterRefKey<Order, Layer>>();

  py::class_<ClusterRefKey<Order, Layer>, ClusterRefBase>
                        (m, cluster_parent_name.c_str());
}

template<size_t Order, size_t Layer, size_t Layer_nd>
struct add_cluster_refs {
  //! starts recursion
  static void static_for(py::module & m) {
    add_cluster_ref<Order, Layer>(m);
    add_cluster_refs<Order, Layer+1, Layer_nd>::static_for(m);
  }
};

//! Stop the recursion
template<size_t Order, size_t Layer_nd>
struct add_cluster_refs<Order, Layer_nd, Layer_nd> {
    static void static_for(py::module &) {}
};

//! templated function for adding clusters of different orders
template<size_t Order, typename SMI>
decltype(auto) add_cluster(py::module & m) {
  using ClusterRef = ClusterRef_t<SMI, Order>;
  using ClusterRefKey = ClusterRefKey_t<SMI, Order>;

  std::string cluster_name =
                internal::GetBindingTypeName<SMI>();

  // TODO(felix) switch case with throw
  if (Order == 1) {
    cluster_name += std::string(".Center");
  } else if (Order == 2) {
    cluster_name += std::string(".Pair");
  } else if (Order == 3) {
    cluster_name += std::string(".Triplet");
  } else if (Order == 4) {
    cluster_name += std::string(".Quadruplet");
  }

  py::class_<ClusterRef, ClusterRefKey>
                      py_cluster(m, cluster_name.c_str());
  py_cluster.def_property_readonly("atom_index", & ClusterRef::get_atom_index,
                                   py::return_value_policy::reference)
    .def_property_readonly("atom_type", & ClusterRef::get_atom_type,
                           py::return_value_policy::reference)
    .def_property_readonly("index", & ClusterRef::get_index,
                           py::return_value_policy::reference)
    .def_property_readonly("size", & ClusterRef::size,
                           py::return_value_policy::reference)
    .def_property_readonly("position", & ClusterRef::get_position,
                           py::return_value_policy::reference);
  return py_cluster;
}

//! bind iterator and ClusterRef for Order >= 2
template<typename StructureManagerImplementation, size_t Order>
decltype(auto) add_iterator(py::module & m,
      PyClusterRef<StructureManagerImplementation, Order-1>& py_cluster) {
  // Order-1 for PyClusterRef because it is the
  // cluster from the previous iteration
  using Child = StructureManagerImplementation;
  using Parent = typename Child::Parent;

  // bind the iteration over clusterRef<Order-1>
  using ClusterRef = typename Parent::template ClusterRef<Order-1>;
  py_cluster.def("__iter__", [] (ClusterRef & v) {
    return py::make_iterator(v.begin(), v.end());
  }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
  auto py_cluster_new = add_cluster<Order, Child>(m);
  return py_cluster_new;
}

//! bind iterator and ClusterRef for Order == 1
template<typename StructureManagerImplementation, size_t Order>
decltype(auto) add_iterator(py::module & m,
                PyManager<StructureManagerImplementation>& manager) {
  using Child = StructureManagerImplementation;
  using Parent = typename Child::Parent;

  // bind the iteration over the manager
  manager.def("__iter__", [] (Parent & v) {
    return py::make_iterator(v.begin(), v.end());
  }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
  auto py_cluster = add_cluster<1, Child>(m);
  return py_cluster;
}

/**
 * Bind the clusterRef allowing to iterate over the manager, atom, neigh...
 * Use signature overloading to dispatch to the proper function.
 * Use iteration by recursion to iterate from Order to To-1 staticaly
*/
template<typename StructureManagerImplementation, size_t Order, size_t MaxOrder>
struct add_iterators {
  //! starts recursion
  static void static_for(py::module & m,
                PyManager<StructureManagerImplementation> & manager) {
    auto py_cluster =
          add_iterator<StructureManagerImplementation, Order>(m, manager);
    add_iterators<
      StructureManagerImplementation, Order+1, MaxOrder>::static_for(m, py_cluster);
  }
  //! following recursion
  static void static_for(py::module & m,
            PyClusterRef<StructureManagerImplementation, Order-1>& py_cluster) {
    auto py_cluster_new =
          add_iterator<StructureManagerImplementation, Order>(m, py_cluster);
    add_iterators<
      StructureManagerImplementation, Order+1, To>::static_for(m,
                                                               py_cluster_new);
  }
};

//! Stop the recursion
template<typename StructureManagerImplementation, size_t To>
struct add_iterators<StructureManagerImplementation, To, To> {
    static void static_for(py::module & ,
                PyClusterRef<StructureManagerImplementation, To-1>& ) {}
};

//! templated function for adding a StructureManager interface
template<typename StructureManagerImplementation>
decltype(auto) add_structure_manager_interface(py::module & m) {
  using Child = StructureManagerImplementation;
  using Parent = typename Child::Parent;

  std::string manager_name = "StructureManager.";
  manager_name += internal::GetBindingTypeName<Child>();
  py::class_<Parent> manager(m, manager_name.c_str());
  manager.def(py::init<>());
  return manager;
}

//! templated function for adding a StructureManager implementation
template<typename StructureManagerImplementation>
decltype(auto) add_structure_manager_implementation(py::module & m,
                                                    py::module & m_garbage) {
  using Child = StructureManagerImplementation;
  using Parent = typename Child::Parent;
  constexpr static size_t MaxOrder = Child::traits::MaxOrder;

  std::string manager_name =
              internal::GetBindingTypeName<Child>();
  py::class_<Child, Parent>  manager(m, manager_name.c_str());
  manager.def(py::init<>());

  // MaxOrder+1 because it stops at Val-1
  add_iterators<Child, 1, MaxOrder+1>::static_for(m_garbage, manager);
  return manager;
}

//! templated function for adding an adaptor
template<typename Adaptor, typename ...ConstructorPack>
decltype(auto) add_adaptor(py::module & m, py::module & m_garbage) {
  using Child = Adaptor;
  using Parent = typename Child::Parent;
  using Implementation_t = typename Child::Implementation_t;
  constexpr static size_t MaxOrder = Child::traits::MaxOrder;

  std::string adaptor_name =
              internal::GetBindingTypeName<Child>();
  py::class_<Child, Parent>  adaptor(m, adaptor_name.c_str());
  // bind constructor making sure the input Manager stays alive
  // as long as the class is alive
  adaptor.def(py::init<Implementation_t&, ConstructorPack...>(),
              py::keep_alive<1, 2>());
  adaptor.def("update", [](Child& v){v.update();});
  // bind clusterRefs so that one can loop over adaptor
  // MaxOrder+1 because recursion stops at Val-1
  add_iterators<Child, 1, MaxOrder+1>::static_for(m_garbage, adaptor);
  return adaptor;
}

//! Template overloading of the binding of the structure managers
template<typename StructureManagerImplementation>
void add_structure_manager(py::module & mod, py::module & m_garbage);

template<typename StructureManagerCenters>
void add_structure_manager(py::module & mod, py::module & m_garbage) {
  using Manager_t = StructureManagerCenters;

  // bind parent class
  add_structure_manager_interface<Manager_t>(m_garbage);
  // bind implementation class
  auto manager =
       add_structure_manager_implementation<Manager_t>(mod, m_garbage);
  manager.def("update", [] (Manager_t & v,
                       const py::EigenDRef<const Eigen::MatrixXd> & positions,
                       const py::EigenDRef<const Eigen::VectorXi> & atom_types,
                       const py::EigenDRef<const Eigen::MatrixXd> & cell,
                       const py::EigenDRef<const Eigen::MatrixXi> & pbc ) {
           v.update(positions, atom_types, cell, pbc);
         });
}

//! Function defining which adaptors are stacked on top of each and then the binding is done for the 
//! structure manager
template<typename StructureManagerImplementation>
void add_adaptors(py::module & mod, py::module & m_garbage) {
  using Manager_t = StructureManagerImplementation;

  using AdaptedManager_0_t = AdaptorNeighbourList<Manager_t>;
  // bind parent class
  add_structure_manager_interface<AdaptedManager_0_t>(m_garbage);
  // bind implementation class
  add_adaptor<AdaptedManager_0_t, double>(mod, m_garbage);

  using AdaptedManager_1_t = AdaptorStrict<AdaptedManager_0_t>;
  // bind parent class
  add_structure_manager_interface<AdaptedManager_1_t>(m_garbage);
  // bind implementation class
  add_adaptor<AdaptedManager_1_t, double>(mod, m_garbage);
}

//! Utility function allowing template differentiation on the structure_manager
template<typename StructureManagerImplementation>
void bind_structure_manager(py::module & m_str_mng, py::module & m_adp,
                            py::module & m_garbage) {
  add_structure_manager<StructureManagerImplementation>(m_str_mng, m_garbage);
  add_adaptors<StructureManagerImplementation>(m_adp, m_garbage);
}

void bind_cluster_ref_base(py::module & m_garbage) {
  py::class_<ClusterRefBase> (m_garbage, "ClusterRefBase");
}

void bind_cluster_refs(py::module & m_garbage) {
  add_cluster_refs<1, 0, 6>::static_for(m_garbage);
  add_cluster_refs<2, 0, 6>::static_for(m_garbage);
  add_cluster_refs<3, 0, 6>::static_for(m_garbage);
  add_cluster_refs<4, 0, 6>::static_for(m_garbage);
}

//! Main function to add StructureManagers and theirs Adaptors
void add_structure_managers(py::module & m_str_mng, py::module & m_adp,
                            py::module & m_garbage) {
  bind_cluster_ref_base(m_garbage);
  bind_cluster_refs(m_garbage);

  bind_structure_manager<StructureManagerCenters>(m_str_mng, m_adp, m_garbage);
}

