#!/usr/bin/env python
"""Example for selecting centres by species

This computes the SphericalInvariants on a small-molecule test set, selecting
only carbon centres
"""

import sys
import ase.io
import numpy as np

sys.path.insert(0,'/local/scratch/source/librascal/build/')
from rascal.representations import SphericalInvariants
from rascal.representations import SphericalInvariantsKspace

#molecules = ase.io.read('/local/scratch/source/librascal/examples/test_lode/water_10angs.xyz',":")
molecules = ase.io.read('/local/scratch/source/librascal/examples/test_lode/water_20angs.xyz',":")
#molecules = ase.io.read('/local/scratch/source/librascal/examples/test_lode/water_50angs.xyz',":")
#molecules = ase.io.read('/local/scratch/source/librascal/examples/test_lode/water_100angs.xyz',":")
hypers = {'interaction_cutoff': 4.0,
          'cutoff_smooth_width': 0.0,
          'max_radial': 8,
          'max_angular': 0,
          'gaussian_sigma_type': "Constant",
          'gaussian_sigma_constant': 0.5}


# direct 
representation = SphericalInvariants(**hypers)
atoms_transformed = representation.transform(molecules)
featvec = atoms_transformed.get_features(representation)
np.savetxt("kernel_direct_20angs.txt",np.dot(featvec,featvec.T))

# reciprocal
representation = SphericalInvariantsKspace(**hypers)
atoms_transformed = representation.transform(molecules)
featvec = atoms_transformed.get_features(representation)
np.savetxt("kernel_reciprocal_20angs.txt",np.dot(featvec,featvec.T))
