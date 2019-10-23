#!/usr/bin/env python
"""Example for selecting centres by species

This computes the SphericalInvariants on a small-molecule test set, selecting
only carbon centres
"""


import ase.io
import numpy as np

from rascal.representations import SphericalInvariants
from rascal.neighbourlist.structure_manager import (
        mask_center_atoms_by_species, mask_center_atoms_by_id)


molecules = ase.io.read('data/small_molecules-1000.xyz', ':100')
n_centers = sum(mol.get_number_of_atoms() for mol in molecules)
n_carbon_centers = sum(np.sum(mol.get_atomic_numbers() == 6)
                       for mol in molecules)
print("{:d} molecules, {:d} centres total, {:d} carbon centres".format(
    len(molecules), n_centers, n_carbon_centers))

for molecule in molecules:
    mask_center_atoms_by_species(molecule, species_select=['C',])
    # Also works by atomic number
    #mask_center_atoms_by_species(molecule, species_select=[6,])

hypers = {'interaction_cutoff': 5.0,
          'cutoff_smooth_width': 0.5,
          'max_radial': 8,
          'max_angular': 6,
          'gaussian_sigma_type': "Constant",
          'gaussian_sigma_constant': 0.3}

representation = SphericalInvariants(**hypers)
atoms_transformed = representation.transform(molecules)
print("Number of feature vectors computed: {:d}".format(
    atoms_transformed.get_dense_feature_matrix(representation).shape[0]))

print("Now masking out the first 5 atoms of each molecule.")
n_remaining_centers = sum(np.sum((mol.get_atomic_numbers()[5:] == 6))
                          for mol in molecules)
print("Number of centres remaining: {:d}".format(n_remaining_centers))

for molecule in molecules:
    mask_center_atoms_by_id(molecule, id_blacklist=np.arange(5))
atoms_transformed = representation.transform(molecules)
print("Number of feature vectors computed: {:d}".format(
    atoms_transformed.get_dense_feature_matrix(representation).shape[0]))
