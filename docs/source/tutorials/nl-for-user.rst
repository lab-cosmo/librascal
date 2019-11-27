.. _nl-for-user:

How to build a neighbor list?
=============================

Overview
~~~~~~~~


Here we provide a tutorial on how to use the C++ interface for building neighbor lists.

Rascal comes with a built-in neighborhood manager. The role of a neighborhood manager is to minimize manual handeling of building the neighborhood lists.
The current implementation avoids the explicit use of embedded lists to describe neighborhoods and makes iterations over them more trivial.
The neighborhood manager interface is universal to all the representations of Rascal. It can also be very easily interfaced with MD/MC simulation packages.
We prefer using `Eigen <http://eigen.tuxfamily.org/>`_ library for the linear algebra operations because of the versatility and flexibility allowed by the library.

`Eigen` library
~~~~~~~~~~~~~~~

`Eigen` is an open source library for linear algebra operatins in C++. It is implemented using expression templates techniques.
We usually use the following templates to describe the positions, the atomic numbers array,..etc

.. code-block:: cpp

    using Vector_t = Eigen::Matrix<double, dim, 1>; // dim = 1, 2 or 3
    Eigen::MatrixXd MATRIX(NROWS,NCOLUMNS); //
    using rascal::VecXi = Eigen::Matrix<int, Eigen::Dynamic, 1>

Creating a neighborhood manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, the neighborhood manager is created and initialized using the :cpp:class:`rascal::NeighborhoodManagerCell`.

.. code-block:: cpp

     using Manager_t = rascal::NeighbourhoodManagerCell;
     Manager_t manager; // initialize the neighborhood manager

While the neighborhood manager can be built for the first time using  :cpp:class:`rascal::NeighborhoodManagerCell::build`,
we recommend updating it using :class:`rascal::NeighborhoodManagerCell::update`. The following inputs should be provided in order:

- the atomic positions matrix column wise as an `Eigen` matrix object
- an array of the atomic numbers of the species
- an array of centers ids
- the cell vectors column wise as an `Eigen` matrix object
- a boolean array of the periodic boundary conditions over the 3 axes
- the maximum radial cutoff (in the units of the cell vectors and atomic positions)

The folllowing code is an example of the inputs that need to be customized according to the user's use, as in to provide
cell vectors, atomic positions and numbers and the ids of the centers.

.. code-block:: cpp

     Eigen::MatrixXd cell(3,3); // the cell vectors matrix

     int number_of_atoms{100}; // the number of atoms in the structure
     Eigen::MatrixXd positions(3,number_of_atoms); // atomic postions matrix
     rascal::VecXi numbers(22); // atomic numbers array

     std::vector<int> center_ids{}; // a vector containing the atomic centers ids

     std::array<bool, 3> pbc{true, false, false}; // for periodic boundary conditions over the x axis only

     double cutoff{3.5}; // a radial cutoff distance of 3.5 angstrom

     manager.update(positions, numbers, center_ids, cell, pbc, cutoff);

- .. note:: A neighborhood manager can be updated as much as necessary by providing all the upmentioned inputs.


Commands of the neighborhood manager:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the current implementation of the neighborhood manager, the centers and their neighbors are :cpp:class:`rascal::Neighborhood Cell::ClusterRef` objects of order 1 and 2, respectively.
In order to access objects of order 2, we need to iterate over the corresponding order 1 object (the corresponding center).
This implementation has the advantage of not using explicit arrays to deal with center and neighbor properties.

Several methods are implemented for the :cpp:class:`rascal::NeighborhoodCell::ClusterRef` objects (centers and their neighbors), such as retrieving
the index of the object, its position and its atom type.

- .. note:: The positions of order 2 are given with a certain offset relative to the position of the corresponding center of order 1. If one atom is included in the neighborhood of two different centers, it will have different positions depending on
            the center being iterated over.

- .. warning:: The relative positions of the neighbors of a center are calculated on the fly and are not stored. If needed they have to be stored manually.

This is an example of code than can be used to print to the screen the positions, types and indices of the centers and their neighbors. It also calculates the relative shift in the positions of every center's neighbors
 and its norm.

.. code-block:: cpp

    for (auto center : manager) {

        int center_index{center.get_atom_index()}; // get the index of the center atomic
        Eigen::MatrixXd position{center.get_position()}; // get the position vector of the center
        auto center_type{center.get_atom_type()}; // get the atome type of the center
        std::cout << "Neighbors properties : " << endl;

        for (auto neigh : center.pairs()){
            int neigh_index{center.get_atom_index()}; // get the index of the neighbor
            Eigen::MatrixXd position{neigh.get_position()}; // get the position vector of the neighbor
            auto neigh_type{center.get_atom_type()}; // get the atome type of the neighbor
            auto relative_shift{position - center.get_position()}; // compute the position offset

            std::cout << "This is the position of atom " << neigh_index << " of a type " << neigh_type << endl;
            std::cout << "The relative position is : " << neigh_position << endl;
            std::cout << "The relative shift is : " << relative_shift << endl;
            std::cout << "The norm of the shif is : " << relative_shift.norm() << " ang" << endl;
        }
    }
