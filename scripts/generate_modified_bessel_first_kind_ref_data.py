import sys,os
import json
import ase
import argparse
from mpmath import mp,besseli,pi,sqrt,exp

import numpy as np
import ubjson
import json

# Computes the sample points and weights for Gauss-Legendre quadrature
# and rescales them.
def get_leggauss(order, a, b):
    x,w = leggauss(order)
    # rescaling
    x = (b-a)*0.5 * x + 0.5*(a+b)
    w = (b-a)*0.5 * w
    return x,w

mp.dps = 20; mp.prec = 100;

def sbesseli(n,z):
    """e^{-x}*i_n(x)"""
    return sqrt(pi/(2*z))*besseli(n+0.5,z)*exp(-z)
def sbesseli_complete_square(n,a,r,x):
    """i_n(2arx)*\exp{-ar^2}*\exp{-ax^2}"""
    z = 2*a*r*x
    return float(sqrt(pi/(2*z))*besseli(n+0.5,z)*exp(-a*r**2)*exp(-a*x**2))


def dump_reference_json():
    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))
    data = dict(i_exp=[],i_complete_square=[])
    max_order = 20
    orders = list(range(max_order))
    xs = np.logspace(-2, 3.8, 300)
    for x in xs:
        vals = []
        for order in orders:
            val = sbesseli(order, x)
            vals.append(float(val))
        data["i_exp"].append(dict(x=x,max_order=max_order,vals=vals))

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
                    val = sbesseli_complete_square(order,alpha,rij,xn)
                    vals[-1].append(val)
            data["i_complete_square"].append(
                    dict(alpha=alpha,rij=rij,xs=xns.tolist(),
                         max_order=max_order,vals=vals))

    # data = [data]
    # with open(
    #       (path +
    #        "tests/reference_data/modified_bessel_first_kind_reference.ubjson",
    #       'wb') as f:
    #     ubjson.dump(data,f)
    with open(os.path.join(
            path, "tests", "reference_data",
            "modified_bessel_first_kind_reference.json"), 'w') as f:
        json.dump(data,f)

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
