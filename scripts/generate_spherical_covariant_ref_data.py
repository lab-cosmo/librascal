import argparse
import ase
from ase.io import read
import json
import sys
import os
import numpy as np

sys.path.insert(0, '../build/')
import rascal
import rascal.lib as lrl
from rascal.utils import ostream_redirect
from rascal.representations import SphericalCovariants

root = os.path.abspath('../')
rascal_reference_path = os.path.join(root, 'reference_data/')
inputs_path = os.path.join(rascal_reference_path, "inputs")
read_inputs_path = os.path.join('reference_data/', "inputs")
dump_path = os.path.join('reference_data/', "tests_only")

#############################################################################


def get_feature_vector(hypers, frames):
    with ostream_redirect():
        soap = SphericalCovariants(**hypers)
        soap_vectors = soap.transform(frames)
        print('Feature vector size: %.3fMB' %
              (soap.get_num_coefficients() * 8.0 / 1.0e6))
        feature_vector = soap_vectors.get_features(soap)
    return feature_vector

#############################################################################


def dump_reference_json():
    import ubjson
    from copy import copy
    from itertools import product
    sys.path.insert(0, os.path.join(root, 'build/'))
    sys.path.insert(0, os.path.join(root, 'tests/'))

    cutoffs = [3]
    gaussian_sigmas = [0.4]
    max_radials = [1, 2, 3]
    max_angulars = [1, 2, 4]
    soap_types = ["LambdaSpectrum"]
    inversion_symmetries = [True, False]

    Lambdas = [1]
    fns = [
        os.path.join(inputs_path, "CaCrP2O7_mvc-11955_symmetrized.json"),
        os.path.join(inputs_path, "small_molecule.json")
    ]
    fns_to_write = [
        os.path.join(read_inputs_path, "CaCrP2O7_mvc-11955_symmetrized.json"),
        os.path.join(read_inputs_path, "small_molecule.json"),
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
                          "covariant_lambda": Lambda}
                soap = SphericalCovariants(**hypers)
                soap_vectors = soap.transform(frames)
                x = soap_vectors.get_features(soap)
                x[np.abs(x) < 1e-300] = 0.
                data['rep_info'][-1].append(dict(feature_matrix=x.tolist(),
                                                 hypers=copy(soap.hypers)))

    with open(os.path.join(root, dump_path, "spherical_covariants_reference.ubjson"),
              'wb') as f:
        ubjson.dump(data, f)

##############################################################################


def main(json_dump, save_kernel):

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
