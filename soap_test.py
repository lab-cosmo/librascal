#!/usr/bin/env python3

# Try to dump SOAP vectors (features) for each atom of each molecule

import sys
import numpy as np
from ase.io import read as ase_read
from rascal.representations import SphericalInvariants as SOAP

# What is this integer? Unit? Reasonable range to explore?
max_radial = 6
# What is this integer? Unit? Reasonable range to explore?
max_angular = 6
# what is a power spectrum?
soap_type = "PowerSpectrum"
#           "RadialSpectrum" # only other possibility
if soap_type == "RadialSpectrum":
    max_angular = 0

input_fn = "reference_data/inputs/small_molecules-20.json"

# only those atomic numbers in input_fn
global_species = [1,6,7,8]

# doc can be found at
# https://github.com/lab-cosmo/librascal/bindings/rascal/representations/
#   spherical_invariants.py

hyper_params = {
    # Angstroms? what is the reasonable range to explore for this value?
    "interaction_cutoff": 3.5,
    # Angstroms? what is the reasonable range to explore for this value?
    "cutoff_smooth_width": 0.5,
    "max_radial": max_radial,
    "max_angular": max_angular,
    "gaussian_sigma_type": "Constant", # only possible choice
    # Angstroms? what is the reasonable range to explore for this value?
    "gaussian_sigma_constant": 0.5,
    "soap_type": soap_type,
    # 'RadialScaling' is the only other choice
    "cutoff_function_type": "ShiftedCosine",
    "normalize": True, # recommended
    # 'DVR' is the only other choice
    "radial_basis": "GTO",
    # # as far as I can tell, global_species has no effect on fixing the
    # # dimensionality of the SOAP feature vectors
    # "global_species": global_species
}
soap = SOAP(**hyper_params)

read_all = ":" # ASE constant meaning "read all molecules"
molecules = ase_read(input_fn, read_all)

# FBR: we cannot put too many molecules at once in memory; or it will start
#      swapping; the problem is that SOAP features are high dimensional

for molecule in molecules:
    m = molecule
    m.cell = [200,200,200] # big enough unit cell
    m.positions[:] += 100 # center atoms in this unit cell
    m.pbc = [False, False, False] # disable periodic boundary conditions
    soap_encoded_molecule = soap.transform(m)
    # n_atoms x n_features
    soap_np_array = soap_encoded_molecule.get_features(soap)
    nb_atoms = soap_np_array.shape[0]
    nb_features = soap_np_array.shape[1]
    print("Data matrix: (%d, %d)" % soap_np_array.shape, file=sys.stderr)
    # dump each atom's index; followed by its SOAP features
    for i_atom in range(nb_atoms):
        print('%d' % i_atom, end='')
        for j_feat in range(nb_features):
            print(' %d:%f' % (j_feat, soap_np_array[i_atom][j_feat]), end='')
        print()
