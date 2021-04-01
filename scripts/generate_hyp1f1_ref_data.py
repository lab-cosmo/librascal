
import sys
import os
import json
import ase
import argparse
from mpmath import mp, hyp1f1
import ubjson
import numpy as np

# dump radial and power spectra for methane

root = os.path.abspath('../')
rascal_reference_path = os.path.join(root, 'reference_data/')
inputs_path = os.path.join(rascal_reference_path, "inputs")
dump_path = os.path.join('reference_data/', "tests_only")


def dump_reference_json():
    sys.path.insert(0, os.path.join(root, 'build/'))
    sys.path.insert(0, os.path.join(root, 'tests/'))
    mp.dps = 200
    data = []
    for l in range(20):
        for n in range(20):
            # z > 660 will lead to larger than double::max() values for a > 19
            for z in [-1e2, -1, -1e-2, 0, 1e-2, 1e-1, 1, 10, 20, 30, 40, 60, 80, 100, 150, 200, 500, 660]:
                a = 0.5 * (n + l + 3)
                b = l + 1.5
                # if b > 0 and (np.abs(a) + 1) * np.abs(z) > 0.9 * b and z < 0:
                #     continue
                # if z < 0 and a-b < 1: continue
                val = float(hyp1f1(a, b, z))
                der = float(a / b * hyp1f1(a + 1, b + 1, z))
                data.append(dict(a=a, b=b, z=z, val=val, der=der))
    # a = 5.
    # b = 1.5
    # z = -4.231593
    # val = float(hyp1f1(a, b, z))
    # der = float(a / b * hyp1f1(a + 1, b + 1, z))
    # data.append(dict(a=a, b=b, z=z, val=val, der=der))
    print(len(data))
    with open(os.path.join(root, dump_path,
                           "hyp1f1_reference.ubjson"), 'wb') as f:
        ubjson.dump(data, f)

##########################################################################################
##########################################################################################


def main(json_dump):
    if json_dump == True:
        dump_reference_json()

##########################################################################################
##########################################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_dump', action='store_true',
                        help='Switch for dumping json')

    args = parser.parse_args()
    main(args.json_dump)
