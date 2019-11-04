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
