
import sys,os
import json
import argparse
import numpy as np
import ubjson
import json
from scipy.special import sph_harm, lpmn

# shape (nb_unit_vectors,3)
def load_unit_vectors_from_json():
    fn = "../tests/reference_data/spherical_harmonics_test.json"
    with open(fn,'r') as f:
        data = json.load(f)
    return data['unit_vectors']

# Sized (l_max+0)**2, contains the l,m components in compressed
# format, i.e. (00)(1-1)(10)(11)(2-2).... as it is in our spherical
# harmonics calculation

def get_ascending_angular_lists(max_angular_l):
    ascending_angular_m_list = []
    ascending_angular_l_list = []
    for angular_l in range(max_angular_l+1):
        ascending_angular_m_list += list(range(-angular_l,angular_l+1))
        ascending_angular_l_list += [angular_l] * (angular_l*2+1)
    return ascending_angular_m_list, ascending_angular_l_list


def dump_lpmn_reference_json():
    # TODO(alex) lpmn think about how the matrix should be saved because they are a (m+1,n+1) array
    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))
    data = []

    directions = load_unit_vectors_from_json()
    max_angular_l = 31
    for direction in directions:
        lpmn_results = []
        z = direction[2] # = cos_theta
        # one could use numpy broadcasting, but this is more readable
        for angular_l in range(max_angular_l+1):
            for angular_m in range(-angular_l, angular_l+1):
                result = np.real(lpmn(angular_m, angular_l, z))
                lpmn_results.append(result)
        data.append(dict(max_angular_l=max_angular_l, direction=direction, associated_legendre=lpmn_results))
    print(len(data))
    with open(path+"tests/reference_data/associated_legendre_reference.ubjson",'wb') as f:
        ubjson.dump(data,f)
    return

def dump_reference_json():
    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))
    data = []

    directions = load_unit_vectors_from_json()
    # double cos_theta = direction[2];
    #theta = np.arccos(directions[:,2])
    # double phi = std::atan2(direction[1], direction[0]);
    #phi = np.arctan2(directions[:,1], directions[:,0])
    max_angular_l = 31
    for direction in directions:
        sph_harm_results = []
        theta = np.arccos(direction[2])
        phi = np.arctan2(direction[1], direction[0])
        # one could use numpy broadcasting, but this is more readable
        for angular_l in range(max_angular_l+1):
            for angular_m in range(-angular_l, angular_l+1):
                result = np.real(sph_harm(angular_m, angular_l, theta, phi))
                sph_harm_results.append(result)
        data.append(dict(max_angular_l=max_angular_l, direction=direction, spherical_harmonics=sph_harm_results))
    print(len(data))
    with open(path+"tests/reference_data/spherical_harmonics_reference.ubjson",'wb') as f:
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
