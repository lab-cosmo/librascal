import argparse
import ase
from ase.io import read
import json
import sys
import numpy as np

sys.path.insert(0, '../build/')
import rascal.lib as lrl
import rascal
from rascal.utils import ostream_redirect
from rascal.representations import SphericalCovariants

rascal_reference_path = 'reference_data/'
inputs_path = os.path.join(rascal_reference_path, "inputs")
dump_path = os.path.join(rascal_reference_path, "tests_only")

#############################################################################


def get_feature_vector(hypers, frames):
    with ostream_redirect():
        soap = SphericalCovariants(**hypers)
        soap_vectors = soap.transform(frames)
        print('Feature vector size: %.3fMB' %
              (soap.get_num_coefficients()*8.0/1.0e6))
        feature_vector = soap_vectors.get_features(soap)
    return feature_vector

#############################################################################


def dump_reference_json():
    import ubjson
    import os
    from copy import copy
    from itertools import product
    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))

    cutoffs = [3]
    gaussian_sigmas = [0.4]
    max_radials = [1, 2, 3]
    max_angulars = [1, 2, 4]
    soap_types = ["LambdaSpectrum"]
    inversion_symmetries = [True, False]

    Lambdas = [1]
    fns = [
        os.path.join(
            path, inputs_path, "CaCrP2O7_mvc-11955_symmetrized.json"),
        os.path.join(path, inputs_path, "small_molecule.json")
    ]
    fns_to_write = [
        os.path.join(dump_path, "CaCrP2O7_mvc-11955_symmetrized.json"),
        os.path.join(dump_path, "small_molecule.json"),
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
            for (soap_type, gaussian_sigma, max_radial,
                 max_angular, inversion_symmetry, Lambda) in product(
                    soap_types, gaussian_sigmas, max_radials, max_angulars,
                    inversion_symmetries, Lambdas):
                hypers = {"interaction_cutoff": cutoff,
                          "cutoff_smooth_width": 0.5,
                          "max_radial": max_radial,
                          "max_angular": max_angular,
                          "gaussian_sigma_type": "Constant",
                          "gaussian_sigma_constant": gaussian_sigma,
                          "cutoff_function_type": "ShiftedCosine",
                          "radial_basis": "GTO",
                          "normalize": True,
                          "soap_type": soap_type,
                          "inversion_symmetry": inversion_symmetry,
                          "lam": Lambda}
                soap = SphericalCovariants(**hypers)
                soap_vectors = soap.transform(frames)
                x = soap_vectors.get_features(soap)
                x[np.abs(x) < 1e-300] = 0.
                data['rep_info'][-1].append(dict(feature_matrix=x.tolist(),
                                                 hypers=copy(soap.hypers)))

    with open(path +
              "tests/reference_data/spherical_covariants_reference.ubjson",
              'wb') as f:
        ubjson.dump(data, f)

##############################################################################


def main(json_dump, save_kernel):

    nmax = 9
    lmax = 5
    lam = 2

    test_hypers = {"interaction_cutoff": 3.0,
                   "cutoff_smooth_width": 0.0,
                   "max_radial": nmax,
                   "max_angular": lmax,
                   "gaussian_sigma_type": "Constant",
                   "gaussian_sigma_constant": 0.3,
                   "soap_type": "LambdaSpectrum",
                   "lam": lam,
                   "inversion_symmetry": True}

    nstr = '2'  # number of structures
    frames = read(os.path.join(inputs_path, 'water_rotations.xyz'), ':'+str(nstr))
    species = set(
        [atom for frame in frames for atom in frame.get_atomic_numbers()])
    nspecies = len(species)
    ncen = np.cumsum([len(frame) for frame in frames])[-1]

    x = get_feature_vector(test_hypers, frames)
    x0 = x.shape[0]
    x = x.reshape((x0, 3, -1, (2*lam + 1), nmax**2))
    x = x.transpose((0, 3, 1, 2, 4))
    x = x.reshape((x0*(2*lam + 1), -1))
    kernel = np.dot(x, x.T)
    kernel = kernel.reshape((x0, (2*lam + 1), x0, (2*lam + 1)))
    kernel = kernel.transpose((0, 2, 1, 3))
    sqrtnorm = np.zeros((x0))
    for i in range(x0):
        sqrtnorm[i] = np.sqrt(np.linalg.norm(kernel[i, i]))
    for i in range(x0):
        for j in range(x0):
            kernel[i, j] /= sqrtnorm[i]*sqrtnorm[j]
    if save_kernel is True:
        np.save(os.path.join(dump_path, 'kernel_soap_example_lambda', str(lam), '.npy'), kernel)

#-------------------dump json reference data------------------------#

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
