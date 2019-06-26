#!/usr/bin/python3

import sys
sys.path.insert(0,'../build/')
import json
import ase
import argparse
import rascal
import rascal.lib as lrl
import numpy as np
from ase.io import read
from rascal.representations import SOAPCovariant
from rascal.utils import ostream_redirect
def load_json(fn):
    with open(fn,'r') as f:
        data = json.load(f)
    return data[str(data['ids'][0])]

def json2ase(f):
    return ase.Atoms(**{v:f[k] for k,v in
dict(positions='positions',atom_types='numbers',pbc='pbc',cell='cell').items()
})


##########################################################################################
##########################################################################################

def get_feature_vector(hypers, frames):
    with ostream_redirect():
        soap = SOAPCovariant(**hypers)
        soap_vectors = soap.transform(frames)
        print('Feature vector size: %.3fMB' % (soap.get_num_coefficients()*8.0/1.0e6))
        feature_vector = soap_vectors.get_feature_matrix()
    return feature_vector

##########################################################################################

def normalise(feature_vector):
    x = feature_vector
    ncen = feature_vector.shape[0]
    for i in range(ncen):
        norm = np.linalg.norm(x[i])
        if norm >= 1.0e-20: x[i] /= norm
    return x

##########################################################################################

#dump radial and power spectra for methane
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
    soap_types = ["LambdaSpectrum"]

    fns = [
        os.path.join(path,"tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json"),
        os.path.join(path,"tests/reference_data/small_molecule.json")
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
        frames = [json2ase(load_json(fn))]
        for cutoff in cutoffs:
            print(fn,cutoff)
            data['rep_info'].append([])
            for soap_type in soap_types:
                for gaussian_sigma in gaussian_sigmas:
                    for max_radial in max_radials:
                        for max_angular in max_angulars:
                            hypers = {"interaction_cutoff": cutoff,
                                      "cutoff_smooth_width": 0.0,
                                      "max_radial": max_radial,
                                      "max_angular": max_angular,
                                      "gaussian_sigma_type": "Constant",
                                      "gaussian_sigma_constant": gaussian_sigma,
                                      "soap_type": soap_type,
                                      "inversion_symmetry" : False,
                                      "lam" : 1}
                            x = get_feature_vector(hypers, frames)
                            data['rep_info'][-1].append(dict(feature_matrix=x.tolist(),
                                                hypers=copy(hypers)))

    with open(path+"tests/reference_data/soap_covariant_reference.ubjson",'wb') as f:
        ubjson.dump(data,f)

##########################################################################################
##########################################################################################

def main(json_dump):

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

    nstr = '2' #number of structures
    frames = read('../tests/reference_data/water_rotations.xyz',':'+str(nstr))
    species = set([atom for frame in frames for atom in frame.get_atomic_numbers()])
    nspecies = len(species)
    ncen = np.cumsum([len(frame) for frame in frames])[-1]

    x = get_feature_vector(test_hypers, frames)
    x = x.T #Eigen column major
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
    np.save('kernel_soap_example_lambda'+str(lam)+'.npy', kernel)

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
