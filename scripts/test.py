import os, sys
from ase.io import read
sys.path.insert(0,"../build/")

import sys
import time
import rascal
import json

import ase
from ase.io import read, write
from ase.build import make_supercell
from ase.visualize import view
import numpy as np
import sys

import json

from rascal.representations import SphericalInvariants as SOAP
from rascal.models import Kernel, PseudoPoints
from rascal.neighbourlist import AtomsList

frames = read('../reference_data/inputs/small_molecules-1000.xyz',':10')
lls = []
for ff in frames:
    lls.append(list(range(len(ff))))

hypers = dict(soap_type="PowerSpectrum",
              interaction_cutoff=3.5,
              max_radial=6,
              max_angular=6,
              gaussian_sigma_constant=0.4,
              gaussian_sigma_type="Constant",
              cutoff_smooth_width=0.5,
              normalize=True,
              radial_basis="DVR",
              compute_gradients=True,
              )
soap = SOAP(**hypers)
zeta=1
kernel = Kernel(soap, name='GAP', zeta=zeta, target_type='Structure', kernel_type='Sparse')

managers = soap.transform(frames)

X_pseudo = PseudoPoints(soap)
X_pseudo.extend(managers, lls)

KNM_down = kernel(managers, X_pseudo, grad=(True, False))

print(KNM_down.shape)