
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
    path = '../' #should be changed
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))

    cutoffs = [2, 3]
    gaussian_sigmas = [0.2, 0.5]
    max_radials = [8, 12]

    fns = [
        os.path.join(path,"tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json"),
        os.path.join(path,"tests/reference_data/small_molecule.json")
    ]
    fns_to_write = [
        "reference_data/small_molecule.json",
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
]

    data = dict(filenames=fns_to_write,
                cutoffs=cutoffs,
                gaussian_sigmas=gaussian_sigmas,
                max_radials=max_radials,
                rep_info=[])

    for fn in fns:
        frames = [json2ase(load_json(fn))]
    for cutoff in cutoffs:
        data['rep_info'].append([])
        for gaussian_sigma in gaussian_sigmas:
            for max_radial in max_radials:
                max_angular = 6

                hypers = {"interaction_cutoff": cutoff,
                          "cutoff_smooth_width": 0.0,
                          "max_radial": max_radial,
                          "max_angular": max_angular,
                          "gaussian_sigma_type": "Constant",
                          "gaussian_sigma_constant": gaussian_sigma}
                x = get_soap_vectors(hypers, frames)
                data['rep_info'][-1].append(dict(feature_matrix=x.tolist(),
                                             hypers=copy(hypers)))

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
