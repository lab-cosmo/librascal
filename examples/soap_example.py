#/usr/bin/python3

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

def get_spectrum(hypers, frames):
    with lrl._rascal.utils.ostream_redirect():
        soap = rascal.representation.SOAP(**hypers)
        soap_vectors = soap.transform(frames)
        spectrum = soap_vectors.get_feature_matrix()
    return spectrum

##########################################################################################

def unravel_power_spectrum(spectrum, nspecies, ncen):
    #rascal exploits the an <-> bn' symmetry without including multiplicty
    #uravel the representation
    spectrum = spectrum.reshape((ncen, int(nspecies*(nspecies+1)/2), -1))
    y = np.zeros(tuple([ncen, nspecies, nspecies] + [i for i in spectrum.shape[2:]]))
    counter = 0
    for j in range(nspecies):
        for k in range(j, nspecies):
            y[:, j, k] = spectrum[:, counter]
            y[:, k, j] = y[:, j, k]
            counter += 1
    y = y.reshape((ncen, -1))
    return y

##########################################################################################

def normalise(spectrum):
    x = spectrum
    ncen = spectrum.shape[0]
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
    path = '/home/willatt/codes/librascal/' #should be changed
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))

    cutoffs = [2, 3]
    gaussian_sigmas = [0.2, 0.3]
    max_radials = [8, 12]
    soap_types = ["RadialSpectrum", "PowerSpectrum"]

    frames = read('../tests/reference_data/methane.xyz',':')
    fns_to_write = ["reference_data/methane.json"]

    data = dict(filenames=fns_to_write,
                cutoffs=cutoffs,
                gaussian_sigmas=gaussian_sigmas,
                max_radials=max_radials,
                soap_types=soap_types)

    #trying to follow the nested list structure of the coulomb matrix reference data
    rep_info = []
#    d = dict(feature_matrices=[],hypers=[])
    for soap_type in soap_types:
        il0 = []
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
                              "gaussian_sigma_constant": gaussian_sigma,
                              "soap_type": soap_type }
                    x = get_spectrum(hypers, frames)
                    d = dict(feature_matrices=[],hypers=[])
                    d['feature_matrices'].append(x.tolist())
                    d['hypers'].append(copy(hypers))
                    il2 += [d]
                il1 += [il2]
            il0 += [il1]
        rep_info += [il0]
    data['rep_info'] = rep_info
    with open(path+"tests/reference_data/soap_reference.ubjson",'wb') as f:
        ubjson.dump(data,f)

##########################################################################################
##########################################################################################

def main(json_dump):

    test_hypers = {"interaction_cutoff": 4.0,
                   "cutoff_smooth_width": 0.0,
                   "max_radial": 25,
                   "max_angular": 6,
                   "gaussian_sigma_type": "Constant",
                   "gaussian_sigma_constant": 0.3,
                   "soap_type": "PowerSpectrum" }

    nmax = test_hypers["max_radial"]
    lmax = test_hypers["max_angular"]
    nstr = '5' #number of structures

    frames = read('../tests/reference_data/dft-smiles_500.xyz',':'+str(nstr))
    species = set([atom for frame in frames for atom in frame.get_atomic_numbers()])
    nspecies = len(species)
    #test_hypers["n_species"] = nspecies #not functional
    ncen = np.cumsum([len(frame) for frame in frames])[-1]

#------------------------------------------nu=1------------------------------------------#

    test_hypers["soap_type"] = "RadialSpectrum"
    x = get_spectrum(test_hypers, frames)
    x = x.T #Eigen column major
    x = normalise(x)
    kernel = np.dot(x, x.T)
    np.save('kernel_soap_example_nu1.npy', kernel)

#------------------------------------------nu=2------------------------------------------#

    test_hypers["soap_type"] = "PowerSpectrum"
    x = get_spectrum(test_hypers, frames)
    x = x.T #Eigen column major
    x = unravel_power_spectrum(x, nspecies, ncen)
    x = normalise(x)
    kernel = np.dot(x, x.T)
    np.save('kernel_soap_example_nu2.npy', kernel)

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
