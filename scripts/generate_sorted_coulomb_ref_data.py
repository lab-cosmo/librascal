"""Script used to generate the sorted_coulomb_reference.ubjson reference file
"""

import os,sys
path = os.path.abspath('/local/git/librascal') # should be changed
sys.path.insert(0, os.path.join(path, 'build/'))
sys.path.insert(0, os.path.join(path, 'tests/'))
import json,ubjson
from copy import copy
from rascal.representation import SortedCoulombMatrix

def load_json(fn):
    with open(fn,'r') as f:
        data = json.load(f)
    return data[str(data['ids'][0])]
def json2ase(f):
    from ase import Atoms
    return Atoms(**{v:f[k] for k,v in
dict(positions='positions',atom_types='numbers',pbc='pbc',cell='cell').items()
})

cutoffs = [2,3,4,5]
sorts = ['row_norm','distance']

fns = [
    os.path.join(path,"tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json"),
    os.path.join(path,"tests/reference_data/small_molecule.json")]
fns_to_write = [
    "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
    "reference_data/small_molecule.json"
]
data = dict(filenames=fns_to_write,cutoffs=cutoffs,rep_info=[])
hypers = dict(central_decay=-1,interaction_cutoff=-1,
              interaction_decay=-1,size=10,sorting_algorithm='')

for fn in fns:
    for cutoff in cutoffs:
        print(fn,cutoff)
        data['rep_info'].append([])
        for sort in sorts:
            rep = SortedCoulombMatrix(cutoff,sorting_algorithm=sort)
            frame = [json2ase(load_json(fn))]
            features = rep.transform(frame)
            test = features.get_feature_matrix()
            hypers['size'] = rep.size
            print(rep.size)
            hypers['sorting_algorithm'] = sort
            data['rep_info'][-1].append(dict(feature_matrix=test.tolist(),
                                             hypers=copy(hypers)))

with open(os.path.join(path,
            "tests/reference_data/sorted_coulomb_reference.ubjson"),'wb') as f:
    ubjson.dump(data,f)