#!/usr/bin/python3.6

#from matplotlib import pyplot as plt
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
nstr = 10

##########################################################################################

frames = ase.io.read('../tests/reference_data/dft-smiles_500.xyz',':'+str(nstr))
species = set([atom for frame in frames for atom in frame.get_atomic_numbers()])
nspecies = len(species)
#test_hypers["n_species"] = nspecies #not functional
ncen = np.cumsum([len(frame) for frame in frames])[-1]

with lrl._rascal.utils.ostream_redirect():
    soap = rascal.representation.SOAP(**test_hypers)
    print("parameters", soap.hypers)

with lrl._rascal.utils.ostream_redirect():
    soap_vectors = soap.transform(frames)
    x = soap_vectors.get_feature_matrix().T

x = x.reshape((ncen, int(nspecies*(nspecies+1)/2), -1))

y = np.zeros(tuple([ncen, nspecies, nspecies] + [i for i in x.shape[2:]]))
counter = 0
for j in range(nspecies):
  for k in range(j, nspecies):
    y[:, j, k] = x[:, counter]
    y[:, k, j] = y[:, j, k]
    counter += 1

y = y.reshape((ncen, -1))

for i in range(ncen):
  norm = np.linalg.norm(y[i])
  if norm > 1.0e-20: y[i] /= norm

#build the kernel
kernel = np.dot(y, y.T)
print(kernel)
np.save('kernel_soap_example.npy', kernel)
