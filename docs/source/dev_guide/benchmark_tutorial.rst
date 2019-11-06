.. _benchmark_tutorial:

.. contents::
   :local:

How to run benchmark for different parameters
---------------------------------------------
To change the parameters of the benchmarks you have to go the corresponding
header file of the benchmark and change the parameters. Here an example for the
interpolator

.. literalinclude:: ../../../performance/benchmarks/benchmark_interpolator.hh
    :language: c++
    :start-after: benchmark-dataset-start
    :end-before: benchmark-dataset-end
    :dedent: 2

At the beginning of the header file there should be a description of each
parameter.

.. literalinclude:: ../../../performance/benchmarks/benchmark_interpolator.hh
    :language: c++
    :start-after: benchmark-parameters-start
    :end-before: benchmark-parameters-end
    :dedent: 2

After changing the parameters, recompile the code and run the benchmark.
Because google benchmark is currently not able to print benchmark results
of type string, string parameters like ``filename`` should be also kept to a
single value to be sure which parameter is used. For other parameters different
combinations can be printed.

How to add a new benchmark 
--------------------------

This file is a skeleton file for writing your own benchmark. It explains the
core featuers of google benchmarks and our own added functionalities in the
``performance/benchmarks/benchmarks.hh`` file.

.. literalinclude:: ../tutorials/benchmark_tutorial.cc
    :language: c++
