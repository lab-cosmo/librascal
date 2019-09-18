

from ase.io import read
import numpy as np
import argparse
import ase
import json
import sys
sys.path.insert(0, '../build/')
import rascal
from rascal.utils import ostream_redirect
from rascal.representations import SphericalInvariants
import rascal.lib as lrl


############################################################################

def get_feature_vector(hypers, frames):
    with ostream_redirect():
        soap = SphericalInvariants(**hypers)
        soap_vectors = soap.transform(frames)
        print('Feature vector size: %.3fMB' %
              (soap.get_num_coefficients()*8.0/1.0e6))
        feature_vector = soap_vectors.get_dense_feature_matrix(soap)
    return feature_vector

##############################################################################


##############################################################################

def dump_reference_json():
    import ubjson
    import os
    from copy import copy
    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))

    cutoffs = [2, 3]
    gaussian_sigmas = [0.2, 0.5]
    max_radials = [4, 10]
    max_angulars = [3, 6]
    soap_types = ["RadialSpectrum", "PowerSpectrum", "BiSpectrum"]
    inversion_symmetry = False

    fns = [
        os.path.join(
            path, "tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json"),
        os.path.join(path, "tests/reference_data/small_molecule.json")
    ]
    fns_to_write = [
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
        "reference_data/small_molecule.json",
    ]

    data = dict(filenames=fns_to_write,
                cutoffs=cutoffs,
                gaussian_sigmas=gaussian_sigmas,
                max_radials=max_radials,
                soap_types=soap_types,
                rep_info=[])

    for fn in fns:
        frames = read(fn)
        for cutoff in cutoffs:
            print(fn, cutoff)
            data['rep_info'].append([])
            for soap_type in soap_types:
                for gaussian_sigma in gaussian_sigmas:
                    for max_radial in max_radials:
                        for max_angular in max_angulars:
                            if 'RadialSpectrum' == soap_type:
                                max_angular = 0
                            if "BiSpectrum" == soap_type:
                                max_radial = 2
                                max_angular = 1
                                inversion_symmetry = True

                            hypers = {"interaction_cutoff": cutoff,
                                    "cutoff_smooth_width": 0.5,
                                    "max_radial": max_radial,
                                    "max_angular": max_angular,
                                    "gaussian_sigma_type": "Constant",
                                    "normalize": True,
                                    "cutoff_function_type": "Cosine",
                                    "radial_basis": "GTO",
                                    "gaussian_sigma_constant": gaussian_sigma,
                                    "soap_type": soap_type,
                                    "inversion_symmetry": inversion_symmetry, }

                            soap = SphericalInvariants(**hypers)
                            soap_vectors = soap.transform(frames)
                            x = soap_vectors.get_dense_feature_matrix(soap)
                            x[np.abs(x) < 1e-300] = 0.
                            data['rep_info'][-1].append(
                                dict(feature_matrix=x.tolist(),
                                     hypers=copy(soap.hypers)))

    with open(path+
                "tests/reference_data/spherical_invariants_reference.ubjson",
                'wb') as f:
        ubjson.dump(data, f)

#############################################################################


def main(json_dump, save_kernel):

    nmax = 9
    lmax = 5
    test_hypers = {"interaction_cutoff": 3.0,
                   "cutoff_smooth_width": 0.0,
                   "max_radial": nmax,
                   "max_angular": lmax,
                   "gaussian_sigma_type": "Constant",
                   "gaussian_sigma_constant": 0.3,
                   "soap_type": "PowerSpectrum"}

    nmax = test_hypers["max_radial"]
    lmax = test_hypers["max_angular"]
    nstr = '2'  # number of structures

    frames = read('../tests/reference_data/dft-smiles_500.xyz', ':'+str(nstr))
    species = set(
        [atom for frame in frames for atom in frame.get_atomic_numbers()])
    nspecies = len(species)
    # test_hypers["n_species"] = nspecies #not functional
    ncen = np.cumsum([len(frame) for frame in frames])[-1]

#--------------------------nu=1------------------------------------------#

    test_hypers["soap_type"] = "RadialSpectrum"
    x = get_feature_vector(test_hypers, frames)
    kernel = np.dot(x, x.T)
    if save_kernel is True:
        np.save('kernel_soap_example_nu1.npy', kernel)

#------------------------------------------nu=2------------------------------#

    test_hypers["soap_type"] = "PowerSpectrum"
    x = get_feature_vector(test_hypers, frames)
    kernel = np.dot(x, x.T)
    if save_kernel is True:
        np.save('kernel_soap_example_nu2.npy', kernel)

#------------------------------------------nu=3-----------------------------#

    frames = read('../tests/reference_data/water_rotations.xyz', ':'+str(nstr))
    species = set(
        [atom for frame in frames for atom in frame.get_atomic_numbers()])
    nspecies = len(species)
    ncen = np.cumsum([len(frame) for frame in frames])[-1]
    nmax = 9
    lmax = 2
    test_hypers["soap_type"] = "BiSpectrum"
    test_hypers["inversion_symmetry"] = False
    test_hypers["max_radial"] = nmax
    test_hypers["max_angular"] = lmax
    x = get_feature_vector(test_hypers, frames)
    kernel = np.dot(x, x.T)
    if save_kernel is True:
        np.save('kernel_soap_example_nu3.npy', kernel)

#------------------dump json reference data--------------------------------#

    if json_dump == True:
        dump_reference_json()

##############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_dump', action='store_true',
                        help='Switch for dumping json')
    parser.add_argument('-save_kernel', action='store_true',
                        help='Switch for dumping json')
    args = parser.parse_args()
    main(args.json_dump, args.save_kernel)
