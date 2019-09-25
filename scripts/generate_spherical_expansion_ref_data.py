import rascal.lib as lrl
import rascal
from rascal.utils import ostream_redirect
from rascal.representations import SphericalExpansion
import json
import ase
import argparse
import numpy as np
from ase.io import read
"""Generate reference data for the librascal spherical expansion"""
import sys
sys.path.insert(0, '../build/')


###############################################################################
###############################################################################

def get_soap_vectors(hypers, frames):
    with ostream_redirect():
        sph_expn = SphericalExpansion(**hypers)
        expansions = sph_expn.transform(frames)
        soap_vectors = expansions.get_dense_feature_matrix(sph_expn)
    return soap_vectors

###############################################################################

# dump spherical expansion


def dump_reference_json():
    import ubjson
    import os
    from copy import copy
    from itertools import product

    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))

    cutoffs = [2, 3]
    gaussian_sigmas = [0.2, 0.5]
    max_radials = [4, 10]
    max_angulars = [3, 6]
    cutoff_smooth_widths = [0., 1.]
    radial_basis = ["GTO", "DVR"]

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
                rep_info=[])

    for fn in fns:
        for cutoff in cutoffs:
            data['rep_info'].append([])
            for (gaussian_sigma, max_radial, max_angular,
                 cutoff_smooth_width, rad_basis) in product(
                    gaussian_sigmas, max_radials, max_angulars,
                    cutoff_smooth_widths, radial_basis):
                frames = read(fn)
                hypers = {"interaction_cutoff": cutoff,
                            "cutoff_smooth_width":
                            cutoff_smooth_width,
                            "max_radial": max_radial,
                            "max_angular": max_angular,
                            "gaussian_sigma_type": "Constant",
                            "cutoff_function_type": "Cosine",
                            "gaussian_sigma_constant":
                            gaussian_sigma,
                            "radial_basis": rad_basis}

                sph_expn = SphericalExpansion(**hypers)
                expansions = sph_expn.transform(frames)
                x = expansions.get_dense_feature_matrix(sph_expn)
                x[np.abs(x) < 1e-300] = 0.
                data['rep_info'][-1].append(
                    dict(feature_matrix=x.tolist(),
                            hypers=copy(sph_expn.hypers)))

    with open(path+"tests/reference_data/spherical_expansion_reference.ubjson",
              'wb') as f:
        ubjson.dump(data, f)

###############################################################################
###############################################################################


def main(json_dump, save_kernel):

    test_hypers = {"interaction_cutoff": 4.0,
                   "cutoff_smooth_width": 0.0,
                   "max_radial": 8,
                   "max_angular": 6,
                   "gaussian_sigma_type": "Constant",
                   "gaussian_sigma_constant": 0.3}

    nmax = test_hypers["max_radial"]
    lmax = test_hypers["max_angular"]
    nstr = '5'  # number of structures

    frames = read('../tests/reference_data/dft-smiles_500.xyz', ':'+str(nstr))
    species = set(
        [atom for frame in frames for atom in frame.get_atomic_numbers()])
    nspecies = len(species)
    ncen = np.cumsum([len(frame) for frame in frames])[-1]

    x = get_soap_vectors(test_hypers, frames)
    if save_kernel is True:
        np.save('spherical_expansion_example.npy', x)

#--------------------------dump json reference data--------------------------#

    if json_dump == True:
        dump_reference_json()

###############################################################################
###############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_dump', action='store_true',
                        help='Switch for dumping json')
    parser.add_argument('-save_kernel', action='store_true',
                        help='Switch for dumping json')
    args = parser.parse_args()
    main(args.json_dump, args.save_kernel)
