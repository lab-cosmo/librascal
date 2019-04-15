#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/property.hh"

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>

using namespace rascal; // NOLINT

/* -------------------- property-typedef-start -------------------- */
// Property<Type, Order, PropertyLayer, NbRow, NbCol>
 
// a scalar value of type size_t assigned to atoms
using AtomScalarProperty = Property<size_t, 1, 0, 1, 1>; 
// a scalar value of type double assigned to pairs of atoms
using PairVectorProperty = Property<double, 2, 0, 1, 1>; 
// a vector with values of type double assigned to pairs of atoms
using AnotherPairVectorProperty = Property<double, 2, 0, 3, 1>; 
// a 2x2 matrix with values of type double assigned to pairs of atoms
using TripleMatrixProperty = Property<double, 3, 0, 2, 2>; 
/* -------------------- property-typedef-end -------------------- */

int main () {
  auto a{1.};
  Eigen::MatrixXd cell(3,3);
  // clang-format off
  cell << a,  0., 0.,
          0., a,  0.,
          0., 0., a;
  // clang-format on

  Eigen::MatrixXd positions(3,4);
  auto p_2 = 0.5 * cell.col(0) + 0.5 * cell.col(1);
  auto p_3 = 0.5 * cell.col(0) + 0.5 * cell.col(2);
  auto p_4 = 0.5 * cell.col(1) + 0.5 * cell.col(2);
  // clang-format off
  positions <<     0.0, p_2[0], p_3[0],
                p_4[0],    0.0, p_2[1],
                p_3[1], p_4[1],    0.0,
                p_2[2], p_3[2], p_4[2];
  // clang-format on
  // same number of atoms, therefore only one necessary


  std::array<int, 3> pbc{{true, true, true}};

  using PBC_t = Eigen::Map<Eigen::Matrix<int, 3, 1>>;

  Eigen::VectorXi atom_types(4);
  atom_types << 1, 1, 1, 1;

  StructureManagerCenters structure_manager{};
  structure_manager.update(positions, atom_types, cell, PBC_t{pbc.data()});
  
  /* -------------------- property-construction-start -------------------- */
  std::string atomic_radius_property_metadata{"size of the atom in pm"};  
  AtomScalarProperty atomic_radius_property{structure_manager,
    atomic_radius_property_metadata};

  // the structure_manager contains only hydrogen atoms, we update
  // the atomic radius porperty by iterating through the structure manager
  size_t hydrogen_radius = 53; // in pm
  // initialise the atomic radius storage
  atomic_radius_property.resize_to_zero();
  for (auto atom : structure_manager) {
    atomic_radius_property.push_back(hydrogen_radius);
  }
  /* -------------------- property-construction-end -------------------- */
  /* -------------------- property-access-start -------------------- */
  for (auto atom : structure_manager) {
    std::cout << "Atom of with atomic number " << atom.get_atom_type()
      << " with atomic radius " << atomic_radius_property[atom] << " pm." 
      << std::endl;
  }
  /* -------------------- property-access-end -------------------- */

  return(0);
}
