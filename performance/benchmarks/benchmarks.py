"""
@file  performance/benchmarks/benchmarks.py

@author  Alexander Goscinski <alexander.goscinski@epfl.ch>

@date   13 October 2019

@brief contains functions for benchmarking the python side

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

import functools
import time

# Global parameters for benchmark
nb_iterations = 5

def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        nb_iterations = kwargs['nb_iterations']
        times = nb_iterations * [0]
        for i in range(nb_iterations):
            start_time = time.perf_counter()
            value = func(*args, **kwargs)
            end_time = time.perf_counter()
            times[i] = end_time - start_time
        run_time = sum(times)/len(times)
        print(f"Finished in {run_time:.6f} secs is mean time for {nb_iterations} runs")
        return value
    return wrapper_timer
