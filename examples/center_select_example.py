#!/usr/bin/env python
"""Example for selecting centres by species

This computes the SphericalInvariants on a small-molecule test set"""

import ase.io
from rascal.representations import SphericalInvariants
from rascal.neighbourlist.structure_manager import add_center_atoms_mask_species

molecules = ase.io.read('data/small_molecules-1000.xyz', ':')

for molecule in molecules:
    add_center_atoms_mask_species(molecule, species_select=['C',])
    # Also works by atomic number
    #add_center_atoms_mask_species(molecule, species_select=[6,])

hypers = {'interaction_cutoff': 5.0,
          'cutoff_smooth_width': 0.5,
          'max_radial': 8,
          'max_angular': 6,
          'gaussian_sigma_type': "Constant",
          'gaussian_sigma_constant': 0.3}

representation = SphericalInvariants(**hypers)
soap_vectors = representation.transform(molecules)
print(soap_vectors.get_dense_feature_matrix().shape)
