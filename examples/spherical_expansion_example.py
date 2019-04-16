#!/usr/bin/python3

import sys
sys.path.insert(0,'../build/')
sys.path.insert(0,'../build/bindings/')
import json
import ase
import argparse
import rascal
import rascal.lib as lrl
import numpy as np
from ase.io import read

##########################################################################################
##########################################################################################

def get_soap_vectors(hypers, frames):
    with lrl._rascal.utils.ostream_redirect():
        sph_expn = rascal.representation.SphericalExpansion(**hypers)
        expansions = sph_expn.transform(frames)
        soap_vectors = expansions.get_feature_matrix()
    return soap_vectors

##########################################################################################

#dump spherical expansion
def dump_reference_json():
    import ubjson
    import os
    from copy import copy
    path = '/home/willatt/codes/librascal/' #should be changed
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))

    cutoffs = [2, 3]
    gaussian_sigmas = [0.2, 0.3]
    max_radials = [8, 12]

    frames = read('../tests/reference_data/methane.xyz',':')
    fns_to_write = ["reference_data/methane.json"]

    data = dict(filenames=fns_to_write,
                cutoffs=cutoffs,
                gaussian_sigmas=gaussian_sigmas,
                max_radials=max_radials)

    #trying to follow the nested list structure of the coulomb matrix reference data
    rep_info = []
    for cutoff in cutoffs:
        il1 = []
        for gaussian_sigma in gaussian_sigmas:
            il2 = []
            for max_radial in max_radials:
                hypers = {"interaction_cutoff": cutoff,
                          "cutoff_smooth_width": 0.0,
                          "max_radial": max_radial,
                          "max_angular": 0,
                          "gaussian_sigma_type": "Constant",
                          "gaussian_sigma_constant": gaussian_sigma}
                x = get_soap_vectors(hypers, frames)
                d = dict(feature_matrices=[],hypers=[])
                d['feature_matrices'].append(x.tolist())
                d['hypers'].append(copy(hypers))
                il2 += [d]
            il1 += [il2]
        rep_info += [il1]
    data['rep_info'] = rep_info
    with open(path+"tests/reference_data/spherical_expansion_reference.ubjson",'wb') as f:
        ubjson.dump(data,f)

##########################################################################################
##########################################################################################

def main(json_dump):

    test_hypers = {"interaction_cutoff": 4.0,
                   "cutoff_smooth_width": 0.0,
                   "max_radial": 8,
                   "max_angular": 6,
                   "gaussian_sigma_type": "Constant",
                   "gaussian_sigma_constant": 0.3}

    nmax = test_hypers["max_radial"]
    lmax = test_hypers["max_angular"]
    nstr = '5' #number of structures

    frames = read('../tests/reference_data/dft-smiles_500.xyz',':'+str(nstr))
    species = set([atom for frame in frames for atom in frame.get_atomic_numbers()])
    nspecies = len(species)
    ncen = np.cumsum([len(frame) for frame in frames])[-1]

    x = get_soap_vectors(test_hypers, frames)
    np.save('spherical_expansion_example.npy', x)

#--------------------------------dump json reference data--------------------------------#

    if json_dump == True:
        dump_reference_json()

##########################################################################################
##########################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_dump', action='store_true', help='Switch for dumping json')
    args = parser.parse_args()
    main(args.json_dump)
