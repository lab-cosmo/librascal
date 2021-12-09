# @file  performance/benchmarks/benchmark_spherical_representations.py
#
# @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
#
# @date   13 October 2019
#
# @brief benchmarks for spherical representations
#
# Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

# This file should be executed from the building directory
# quick solution for loading
import sys
from ase.io import read

import benchmarks

try:
    from rascal.representations import SphericalInvariants
    from rascal.representations import SphericalExpansion
except ImportError:
    sys.exit("Could not import rascal.")


# Benchmark hyperparameters
NB_ITERATIONS_PER_REPRESENTATION = 20
# global variable to make sure that range of interpolator agrees with cutoff
INTERACTION_CUTOFF = 4

input_files = [
    "reference_data/inputs/small_molecule.json",
    "reference_data/inputs/methane.json",
]

optimizations_args = [
    {"type": "None"},
    {
        "type": "Spline",
        "accuracy": 1e-8,
        "range": (0.0001, INTERACTION_CUTOFF),
    },
    {
        "type": "Spline",
        "accuracy": 1e-10,
        "range": (0.0001, INTERACTION_CUTOFF),
    },
    {
        "type": "Spline",
        "accuracy": 1e-12,
        "range": (0.0001, INTERACTION_CUTOFF),
    },
]

radial_bases = ["GTO", "DVR"]


@benchmarks.bench
def transform_representation(representation, frames, **kwargs):
    representation.transform(frames)


def benchmark_spherical_representations(frames, optimization_args, radial_basis):
    hypers = {
        "interaction_cutoff": INTERACTION_CUTOFF,
        "max_radial": 8,
        "max_angular": 6,
        "gaussian_sigma_constant": 0.5,
        "gaussian_sigma_type": "Constant",
        "cutoff_smooth_width": 0.5,
        "radial_basis": radial_basis,
        "optimization": optimization_args,
    }

    print("Timing SphericalExpansion")
    transform_representation(
        SphericalExpansion(**hypers),
        frames,
        nb_iterations=NB_ITERATIONS_PER_REPRESENTATION,
    )

    hypers = {
        "soap_type": "PowerSpectrum",
        "interaction_cutoff": INTERACTION_CUTOFF,
        "max_radial": 8,
        "max_angular": 6,
        "gaussian_sigma_constant": 0.5,
        "gaussian_sigma_type": "Constant",
        "cutoff_smooth_width": 0.5,
        "normalize": False,
        "radial_basis": radial_basis,
        "optimization": optimization_args,
    }
    print("Timing SphericalInvariants")
    transform_representation(
        SphericalInvariants(**hypers),
        frames,
        nb_iterations=NB_ITERATIONS_PER_REPRESENTATION,
    )


if __name__ == "__main__":
    for input_file in input_files:
        frames = read(input_file, index=":")
        for i in range(len(frames)):
            frames[i].wrap(eps=1e-11)
        for radial_basis in radial_bases:
            for optimization_args in optimizations_args:
                print(
                    "Benchmark on file " + input_file + " for optimization",
                    optimization_args,
                    "with basis " + radial_basis,
                )
                benchmark_spherical_representations(
                    frames, optimization_args, radial_basis
                )
                print()
