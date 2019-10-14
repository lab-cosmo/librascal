"""
file  benchmark_spherical_representations.py

@author  Alexander Goscinski <alexander.goscinski@epfl.ch>

@date   13 October 2019

@brief benchmarks for spherical representations

Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)

rascal is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3, or (at
your option) any later version.

rascal is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this software; see the file LICENSE. If not, write to the
Free Software Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

import sys

from ase.io import read

import benchmarks

import rascal
from rascal.representations import SphericalInvariants
from rascal.representations import SphericalExpansion

# This file should be executed from the building directory
# quick solution for loading 
RASCAL_DIR = ".."
BUILD_DIR = RASCAL_DIR+"/build"
sys.path.insert(0, BUILD_DIR)


# Benchmark hyperparameters 
NB_ITERATIONS_PER_REPRESENTATION = 5
input_files = [RASCAL_DIR +'/examples/data/small_molecules-1000.xyz']
optimizations_args = [{'type':'Nothing'}, \
                      {'type':'Spline', 'accuracy':1e-8, \
                       'range': (0,4)}]
radial_bases = ["GTO", "DVR"]

@benchmarks.timer
def transform_representation(representation, frames, **kwargs):
    representation.transform(frames)

def benchmark_spherical_representations(frames, optimization_args, radial_basis):
    hypers = dict(
                  interaction_cutoff=4,
                  max_radial=8,
                  max_angular=6,
                  gaussian_sigma_constant=0.5,
                  gaussian_sigma_type="Constant",
                  cutoff_smooth_width=0.5,
                  radial_basis="GTO",
                  optimization_args=optimization_args
                  )
    print("Timing SphericalExpansion")
    transform_representation(SphericalExpansion(**hypers), frames, nb_iterations=NB_ITERATIONS_PER_REPRESENTATION)

    hypers = dict(soap_type="PowerSpectrum",
                  interaction_cutoff=4,
                  max_radial=8,
                  max_angular=6,
                  gaussian_sigma_constant=0.5,
                  gaussian_sigma_type="Constant",
                  cutoff_smooth_width=0.5,
                  normalize=False,
                  radial_basis="GTO",
                  optimization_args=optimization_args
                  )
    print("Timing SphericalInvariants")
    transform_representation(SphericalInvariants(**hypers), frames, nb_iterations=NB_ITERATIONS_PER_REPRESENTATION)

for input_file in input_files:
    frames = read(input_file,index=":")
    for i in range(len(frames)):
        frames[i].wrap(eps=1e-11)
    for radial_basis in radial_bases:
        for optimization_args in optimizations_args:
            print("Benchmark on file "+input_file+" for optimization_args", optimization_args, "with basis " +radial_basis)
            benchmark_spherical_representations(frames, optimization_args, radial_basis)
