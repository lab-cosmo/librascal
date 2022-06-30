# @file  performance/benchmarks/benchmarks.py
#
# @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
#
# @date   13 October 2019
#
# @brief contains functions for benchmarking the python side
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

import functools
from timeit import default_timer as timer

DEFAULT_ITERATIONS = 100


def bench(func):
    """Print the runtime of the decorated function"""

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        nb_iterations = kwargs.get("nb_iterations", DEFAULT_ITERATIONS)

        total_time = 0.0
        for i in range(nb_iterations):
            start_time = timer()
            value = func(*args, **kwargs)
            end_time = timer()
            total_time += end_time - start_time

        total_time /= nb_iterations
        print(f"Finished in {total_time*1e6:.2f} ns over {nb_iterations} runs")
        return value

    return wrapper_timer
