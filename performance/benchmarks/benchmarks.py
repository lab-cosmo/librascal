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
            start_time = time.perf_counter()    # 1
            value = func(*args, **kwargs)
            end_time = time.perf_counter()      # 2
            times[i] = end_time - start_time    # 3
        run_time = sum(times)/len(times)
        print(f"Finished {func.__name__!r} in {run_time:.6f} secs from {nb_iterations} runs")
        return value
    return wrapper_timer
