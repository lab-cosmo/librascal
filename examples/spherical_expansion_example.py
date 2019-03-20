
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
               "max_radial": 8,
               "max_angular": 6,
               "gaussian_sigma_type": "Constant",
               "gaussian_sigma_constant": 0.3}

nmax = test_hypers["max_radial"]
lmax = test_hypers["max_angular"]
nstr = 100

##########################################################################################

frames = ase.io.read('../tests/reference_data/dft-smiles_500.xyz',':'+str(nstr))
species = set([atom for frame in frames for atom in frame.get_atomic_numbers()])
nspecies = len(species)
ncen = np.cumsum([len(frame) for frame in frames])[-1]

with lrl._rascal.utils.ostream_redirect():
    sph_expn = rascal.representation.SphericalExpansion(**test_hypers)
    print("parameters", sph_expn.hypers)

with lrl._rascal.utils.ostream_redirect():
    expansions = sph_expn.transform(frames)

x = expansions.get_feature_matrix().T

#unravel the representation
coeffs = np.zeros((ncen,nspecies,nmax,lmax+1,2*lmax+1),float)
for icen in range(ncen):
    ii=0
    for ispecies in range(nspecies):
        for l in range(lmax+1):
            for im in range(2*l+1):
                for n in range(nmax):
                    coeffs[icen,ispecies,n,l,im] = x[icen,ii]
                    ii+=1

#build the power spectrum
size = nspecies**2*nmax*nmax*(lmax+1)
power = np.zeros((ncen,size),float)
for icen in range(ncen):
    jj = 0
    for ispecies in range(nspecies):
        for jspecies in range(nspecies):
            for n1 in range(nmax):
                for n2 in range(nmax):
                    for l in range(lmax+1):
                        power[icen,jj] = np.dot(coeffs[icen,ispecies,n1,l],
                                                coeffs[icen,jspecies,n2,l])/ \
                                                np.sqrt(2*l+1)
                        jj+=1
    norm = np.linalg.norm(power[icen])
    if norm >= 1.0e-40: power[icen] /= norm

#build the kernel
kernel = np.dot(power,power.T)
print(kernel)
np.save("kernel_spherical_expansion_example.npy", kernel)
