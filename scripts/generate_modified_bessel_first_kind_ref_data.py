import sys
import os
import json
import ase
import argparse
from mpmath import mp, besseli, pi, sqrt, exp

import numpy as np
import ubjson
import json

root = os.path.abspath('../')
rascal_reference_path = os.path.join(root, 'reference_data/')
inputs_path = os.path.join(rascal_reference_path, "inputs")
dump_path = os.path.join('reference_data/', "tests_only")

mp.dps = 20
mp.prec = 100


def sbesseli(n, z):
    """e^{-x}*i_n(x)"""
    return sqrt(pi/(2*z))*besseli(n+0.5, z)*exp(-z)


def sbesseli_complete_square(n, a, r, x):
    """i_n(2arx)*\exp{-ar^2}*\exp{-ax^2}"""
    z = 2*a*r*x
    return float(sqrt(pi/(2*z))*besseli(n+0.5, z)*exp(-a*r**2)*exp(-a*x**2))


def dump_reference_json():
    sys.path.insert(0, os.path.join(root, 'build/'))
    sys.path.insert(0, os.path.join(root, 'tests/'))
    data = dict(i_complete_square=[])
    # 1 test the special case in the bessel function and 20 test that we
    # can ramp up l_max safely
    max_orders = [1, 20]
    for max_order in max_orders:
        orders = list(range(max_order))
        # gaussian sigma in [0.1, 0.9]
        alphas = np.linspace(0.6, 50, 10)
        # looks at rc up to 10
        xns = np.linspace(0.005, 10, 15)
        # atoms should not be much closer than this
        rijs = np.linspace(0.5, 10, 15)

        for alpha in alphas:
            for rij in rijs:
                vals = []
                for xn in xns:
                    vals.append([])
                    for order in orders:
                        val = sbesseli_complete_square(order, alpha, rij, xn)
                        vals[-1].append(val)
                # avoid values that are too small for ubjson to be interpreted
                # as doubles
                vals = np.array(vals)
                vals[vals < 1e-300] = 0.
                data["i_complete_square"].append(
                    dict(alpha=alpha, rij=rij, xs=xns.tolist(),
                         max_order=max_order, vals=vals.tolist()))

    with open(os.path.join(root, dump_path,
                           "modified_bessel_first_kind_reference.ubjson"), 'wb') as f:
        ubjson.dump(data, f)

###############################################################################
###############################################################################


def main(json_dump):
    if json_dump == True:
        dump_reference_json()

###############################################################################
###############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_dump', action='store_true',
                        help='Switch for dumping json')
    args = parser.parse_args()
    main(args.json_dump)
