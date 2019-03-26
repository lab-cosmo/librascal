#!/usr/bin/python3.6 -u

import sys
sys.path.insert(0,'../build/')
sys.path.insert(0,'../build/bindings/')
import rascal
import rascal.lib as lrl
import json
import ase
from ase.visualize import view
import numpy as np
import sys

##########################################################################################

test_hypers = {"interaction_cutoff": 4.0,
               "cutoff_smooth_width": 0.0, 
               "max_radial": 25,
               "max_angular": 6,
               "gaussian_sigma_type": "Constant",
               "gaussian_sigma_constant": 0.3,
               "soap_type": "PowerSpectrum" }

nmax = test_hypers["max_radial"]
lmax = test_hypers["max_angular"]
nstr = '' #number of structures

##########################################################################################

frames = ase.io.read('../tests/reference_data/dft-smiles_500.xyz',':'+str(nstr))
species = set([atom for frame in frames for atom in frame.get_atomic_numbers()])
nspecies = len(species)
#test_hypers["n_species"] = nspecies #not functional
ncen = np.cumsum([len(frame) for frame in frames])[-1]

#------------------------------------------nu=2------------------------------------------#

with lrl._rascal.utils.ostream_redirect():
    soap = rascal.representation.SOAP(**test_hypers)
    print("parameters", soap.hypers)

with lrl._rascal.utils.ostream_redirect():
    soap_vectors = soap.transform(frames)
    x = soap_vectors.get_feature_matrix().T

x = x.reshape((ncen, int(nspecies*(nspecies+1)/2), -1))

#rascal exploits the an <-> bn' symmetry without including multiplicty
#uravel the representation
y = np.zeros(tuple([ncen, nspecies, nspecies] + [i for i in x.shape[2:]]))
counter = 0
for j in range(nspecies):
    for k in range(j, nspecies):
        y[:, j, k] = x[:, counter]
        y[:, k, j] = y[:, j, k]
        counter += 1

y = y.reshape((ncen, -1))

#normalise the representation
for i in range(ncen):
    norm = np.linalg.norm(y[i])
    if norm > 1.0e-20: y[i] /= norm

#------------------------------------------nu=2------------------------------------------#

#------------------------------------------nu=1------------------------------------------#

test_hypers["soap_type"] = "RadialSpectrum" 
with lrl._rascal.utils.ostream_redirect():
    soap = rascal.representation.SOAP(**test_hypers)
    print("parameters", soap.hypers)

with lrl._rascal.utils.ostream_redirect():
    soap_vectors = soap.transform(frames)
    x = soap_vectors.get_feature_matrix().T

#normalise the representation
y = np.zeros(x.shape)
for i in range(ncen):
    norm = np.linalg.norm(x[i])
    if norm > 1.0e-20: y[i] = x[i]/norm

#------------------------------------------nu=1------------------------------------------#

#build the kernel
kernel = np.dot(y, y.T)
print(kernel)
np.save('kernel_soap_example.npy', kernel)

#----------------------------------------------------------------------------------------#
#--------------------------------dump json reference data--------------------------------#
#----------------------------------------------------------------------------------------#

import ubjson
import os
from copy import copy
path = '/home/willatt/codes/librascal/' #should be changed
sys.path.insert(0, os.path.join(path, 'build/'))
sys.path.insert(0, os.path.join(path, 'tests/'))

cutoffs = [2, 3]
gaussian_sigmas = [0.2, 0.3]
max_radials = [8, 15]

frames = ase.io.read('../tests/reference_data/methane.xyz',':')
fns_to_write = ["reference_data/methane.json"]

data = dict(filenames=fns_to_write,
            cutoffs=cutoffs,
            gaussian_sigmas=gaussian_sigmas,
            max_radials=max_radials)

#trying to follow the nested list structure of the coulomb matrix reference data
rep_info = []
d = dict(feature_matrices=[],hypers=[])
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
                      "soap_type": "RadialSpectrum" }
            with lrl._rascal.utils.ostream_redirect():
                soap = rascal.representation.SOAP(**test_hypers)
                soap_vectors = soap.transform(frames)
                x = soap_vectors.get_feature_matrix()
            d['feature_matrices'].append(x.tolist())
            d['hypers'].append(copy(hypers))
            il2 += [d]
        il1 += il2
    rep_info += il1
data['rep_info'] = rep_info

with open(path+"tests/reference_data/soap_radial_spectrum_reference.ubjson",'wb') as f: 
    ubjson.dump(data,f)

#----------------------------------------------------------------------------------------#
#--------------------------------dump json reference data--------------------------------#
#----------------------------------------------------------------------------------------#
