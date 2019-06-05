import sys
sys.path.insert(0,'../build/')
import json
import ase
import rascal
import numpy as np
from ase.io import read
from rascal.representations import SOAP

N_ITERATIONS = 500

def load_json(fn):
    with open(fn,'r') as f:
        data = json.load(f)
    return data[str(data['ids'][0])]

def json2ase(f):
    return ase.Atoms(**{v:f[k] for k,v in
dict(positions='positions',atom_types='numbers',pbc='pbc',cell='cell').items()
})

fn = '../tests/reference_data/molecular_crystal.json'

frames = [json2ase(load_json(fn))] * N_ITERATIONS

hypers = dict(soap_type="PowerSpectrum",
              interaction_cutoff=5.,
              max_radial=8,
              max_angular=6,
              gaussian_sigma_constant=0.4,
              gaussian_sigma_type="Constant",
              cutoff_smooth_width=0.5,
              )

soap = SOAP(**hypers)


soap.transform(frames)