#!/usr/bin/env python3

# Try to dump SOAP vectors for each atom of a molecule

from ase.io import read as ase_read
from rascal.representations import SphericalInvariants as SOAP

max_radial = 6 # What is this integer? Unit? Reasonable range to explore?
max_angular = 6 # What is this integer? Unit? Reasonable range to explore?
soap_type = "PowerSpectrum"  # only physicists know what is a power spectrum...
#            "RadialSpectrum" # the other possibility; also only known to physicists...
if soap_type == "RadialSpectrum":
    max_angular = 0

input_fn = "reference_data/inputs/small_molecules-20.json"

start = 0
length = 1
representation_name = "spherical_invariants" # Any other choice?
target_type = "Atom" # The other choice is Structure; i.e. granularity is per atom or per molecule

hyper_params = {
    "interaction_cutoff": 3.5, # Angstroms? what is the reasonable range to explore for this value?
    "cutoff_smooth_width": 0.5, # Angstroms? what is the reasonable range to explore for this value?
    "max_radial": max_radial,
    "max_angular": max_angular,
    "gaussian_sigma_type": "Constant", # other possible choices for this value?
    "gaussian_sigma_constant": 0.5, # Angstroms? what is the reasonable range to explore for this value?
    "soap_type": soap_type,
    "cutoff_function_type": "ShiftedCosine", # other possible choices for this value?
    "normalize": True,
    "radial_basis": "GTO", # other possible choices for this value?
}
soap = SOAP(**hyper_params)

molecules = ase_read(input_fn, "{}:{}".format(start, start + length))
soap_encoded_molecules = soap.transform(molecules)

for molecule in soap_encoded_molecules:
    for soap_encoded_atom in molecule:
        print(soap_encoded_atom)
