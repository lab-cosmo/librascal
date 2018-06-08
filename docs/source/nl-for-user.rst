.. _nl-for-user:

How to build a neighbor list?
=================================

Here we provide a tutorial on how to use the C++ interface for building neighbor lists.

Rascal comes with a neighborhood manager. We prefer using `Eigen` library for treating matrices because of the versatility and flexibility allowed by the library.

First, the neighborhood manager is initialized using the :class:`rascal::NeighborhoodManagerCell`.

.. code-block:: cpp

     using Manager_t = rascal::NeighbourhoodManagerCell;
     Manager_t manager; // initialize the neighborhood manager

Then, the  neighborhood manager needs to be updated using :class:`rascal::NeighborhoodManagerCell::update`. The following inputs should be provided: 

- the cell vectors column wise
- the atomic positions matrix column wise
- centers ids
- the periodic boundary conditions over the 3 axes
- the maximum radial cutoff

.. code-block:: cpp

     Eigen::MatrixXd cell(3,3); // the cell vectors matrix
     
     int number_of_atoms{100}; // the number of atoms in the structure
     Eigen::MatrixXd positions(3,number_of_atoms); // atomic postions matrix
     
     std::vector<int> center_ids{}; // a vector containing the atomic centers ids
     
     std::array<bool, 3> pbc{true, false, false}; // for periodic boundary conditions over the x axis only
     
     double cutoff{3.5}; // a radial cutoff distance of 3.5 angstrom

     manager.update(positions, center_ids, cell, pbc, cutoff);


Methods implemented for neighborhood manager:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In the current implementation of the neighborhood manager, the centers and their neighbors are :class:`rascal::Neighborhood Cell::ClusterRef` objects of level 1 and 2, respectively. 
In order to access objects of level 2, we need to iterate over the corresponding leve 1 object.
This implementation has the advantage of not using explicit arrays to deal with center and neighbor properties. 


- .. warning:: The neighbors' properties are only accessible from their centers objects. The relative positions of the neighbors of a center are calculated on the fly and are not stored. 

Several methods are implemented for the :class:`rascal::Neighborhood Cell::ClusterRef` objects (centers and their neighbors), such as retrieving 
the index of the object, its position and its atom type.

- .. note:: The positions of level 2 are given with a certain offset relative to the position of the corresponding center of level 1.

This is an example of code than can be used to print to the screen the positions, types and indices of the centers and their neighbors. It shows also the relative shift in the positions of every center's neighbors.

.. code-block:: cpp

    for (auto center : manager) {
        int center_index{center.get_atom_index()}; // get the index of the center atomic
        Eigen::MatrixXd position{center.get_position()}; // get the position vector of the center
        auto center_type{center.get_atom_type()}; // get the atome type of the center
        
        for (auto neigh : center){
            int neigh_index{center.get_atom_index()}; // get the index of the neighbor 
            Eigen::MatrixXd position{neigh.get_position()}; // get the position vector of the neighbor
            auto neigh_type{center.get_atom_type()}; // get the atome type of the neighbor
            auto relative_shift{position - center.get_position()}; // compute the position offset

            std::cout << "This is the position of atom " << neigh_index << " of a type " << neigh_type << endl;
            std::cout << "The relative position is : " << neigh_position << endl;
            std::cout << "The relative shift is : " << relative_shift << endl;
        }
    }
