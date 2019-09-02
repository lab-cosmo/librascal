
import sys
sys.path.insert(0,'../build/')
import json
import ase
import argparse
import rascal
import rascal.lib as lrl
import numpy as np
from ase.io import read
from rascal.representations import SphericalInvariants
from rascal.models import Kernel
from rascal.utils import ostream_redirect



#dump radial and power spectra for methane
def dump_reference_json():
    import ubjson
    import os
    from copy import copy
    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))

    cutoffs = [3.5]
    gaussian_sigmas = [0.5]
    max_radials = [6]
    max_angulars = [6]
    soap_types = ["RadialSpectrum", "PowerSpectrum"]

    fn = os.path.join(path,"tests/reference_data/dft-smiles_500.xyz")
    fn_to_write = "reference_data/dft-smiles_500.ubjson"
    start = 0
    length = 5
    representations = ['spherical_invariants']
    kernel_names = ['Cosine']
    target_types = ['Structure', 'Atom']
    dependant_args = dict(Cosine=[dict(zeta=1),dict(zeta=2),dict(zeta=4)])

    data = dict(filename=fn_to_write,
                start=start,
                length=length,
                cutoffs=cutoffs,
                gaussian_sigmas=gaussian_sigmas,
                max_radials=max_radials,
                soap_types=soap_types,
                kernel_names=kernel_names,
                target_types=target_types,
                dependant_args=dependant_args,
                rep_info=dict(spherical_invariants=[]))

    frames = read(fn, '{}:{}'.format(start,start+length))
    for representation_name in representations:
        for cutoff in cutoffs:
            print(fn,cutoff)
            data['rep_info'][representation_name].append([])
            for kernel_name in kernel_names:
                for target_type in target_types:
                    for kwargs in dependant_args[kernel_name]:
                        for soap_type in soap_types:
                            for gaussian_sigma in gaussian_sigmas:
                                for max_radial in max_radials:
                                    for max_angular in max_angulars:
                                        if 'RadialSpectrum' == soap_type:
                                            max_angular = 0

                                        hypers = {"interaction_cutoff": cutoff,
                                                "cutoff_smooth_width": 0.5,
                                                "max_radial": max_radial,
                                                "max_angular": max_angular,
                                                "gaussian_sigma_type": "Constant",
                                                "gaussian_sigma_constant": gaussian_sigma,
                                                "soap_type": soap_type,
                                                "cutoff_function_type":"Cosine",
                                                "normalize": True,
                                                "radial_basis":"GTO"}
                                        soap = SphericalInvariants(**hypers)
                                        soap_vectors = soap.transform(frames)
                                        hypers_kernel = dict(name=kernel_name,
                                        target_type=target_type)
                                        hypers_kernel.update(**kwargs)
                                        kernel = Kernel(soap, **hypers_kernel)
                                        kk = kernel(soap_vectors)
                                        # x = get_spectrum(hypers, frames)
                                        for aa in soap.nl_options:
                                            aa['initialization_arguments'] = aa['args']

                                        data['rep_info'][representation_name][-1].append(dict(kernel_matrix=kk.tolist(),
                                        hypers_rep=copy(soap.hypers),
                                        hypers_manager=copy(soap.nl_options),
                                        hypers_kernel=copy(hypers_kernel)))

    with open(path+"tests/reference_data/kernel_reference.ubjson",'wb') as f:
        ubjson.dump(data,f)

##########################################################################################
##########################################################################################

def main(json_dump):


    if json_dump == True:
        dump_reference_json()

##########################################################################################
##########################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_dump', action='store_true', help='Switch for dumping json')
    args = parser.parse_args()
    main(args.json_dump)
