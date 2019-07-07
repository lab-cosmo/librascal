import ase
import os, sys
import numpy as np
sys.path.insert(0,"../build/")

from rascal.representations import SOAPInvariant


def load_json(fn):
    import json
    with open(fn,'r') as f:
        data = json.load(f)
    return data[str(data['ids'][0])]

def json2ase(f):
    return ase.Atoms(**{v:f[k] for k,v in
dict(positions='positions',atom_types='numbers',pbc='pbc',cell='cell').items()
})


fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'

hypers = {'interaction_cutoff': 3,
  'cutoff_smooth_width': 0.0,
  'max_radial': 2,
  'max_angular': 1,
  'gaussian_sigma_type': 'Constant',
  'gaussian_sigma_constant': 0.5,
  'soap_type': 'BiSpectrum',
  'inversion_symmetry': True,
  'normalize': False,
  'radial_basis': 'GTO'}

frames = [json2ase(load_json(fn))]

soap = SOAPInvariant(**hypers)
representation = soap.transform(frames)
X = representation.get_feature_matrix()
print(np.where(np.isnan(X) == True))
print(X[:,0])