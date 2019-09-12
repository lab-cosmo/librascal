.. _how_to_add_representation:

.. contents::
   :local:

How to add a new representation
-------------------------------



Write a Calculator for a representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A calculator of a representation is an object that builds a representation of the atomic structure contained in *a* structure manager. This class:

- inherits publicly from :cpp:class:`CalculatorBase <rascal::CalculatorBase>` to follow its interface and use some of the common utilies shared by such class.

- uses a :cpp:func:`compute() <compute()>` function expecting one or a list of structure manager(s) to build the representation efficiently.

- uses the :cpp:class:`Property <rascal::Property>` class to store the representation's features.

The behaviour of a representation calculator is defined at construction by hyperparameters which are specific to the representation.
. These hyperparameters can
  - define a particular implementation from a set of conceptually equivalent methods, implemented using primitive enums
  - be scalar values or on/off features used as hyperparameters in the calculation of the represenation, implemented using primitive types

Note that there is one representation manager per structure manager.

To illustrate the basic structure that a new representation that would be implemented in ``calculator_representation_name.hh`` should follow, let's take the example of the :cpp:class:`CalculatorSortedCoulomb <rascal::CalculatorSortedCoulomb>`. A detailed discussion of the sorted coulomb matrix representation can be found in [#one]_ and [#two]_.


The representation starts with the definition of some useful short hands

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-preamble-start
    :end-before: rep-preamble-end
    :dedent: 2

followed by the definition of its constructors and destructor

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-construc-start
    :end-before: rep-construc-end
    :dedent: 2

the declaration of the concrete implementation of the calculator interface

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-interface-start
    :end-before: rep-interface-end
    :dedent: 2

and the declaration of some functions for internal use in the protected section. The end of the class contains the different internal variables needed by the class

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-variables-start
    :end-before: rep-variables-end
    :dedent: 2

Implementing several types of the same representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Often there are different ways to implement the same fundamental representation. To be able to do this in rascal, the compute function is used with the type of of the implementation as template argument (or multiple if representation requires this). The switch between several implementations of conceptually equivalent parts of a representation can be implemented through several mechanisms such a virtual inheritance. We detail here how to implement such switch efficiently using the :cpp:class:`CalculatorSortedCoulomb <rascal::CalculatorSortedCoulomb>` as an example.

To make the coulomb matrix invariant with respect to permutation Ref. [#two]_ proposes to sort the upper triangular part of the coulomb matrix according to the norm of each rows or the distance from the central atom (see [#one]_ for details).

The implementation of these two behaviour is encapsulated in the :cpp:class:`SortCoulomMatrix <rascal::internal::SortCoulomMatrix>` class and the choice between them is done with a template parameter using template specialization. Note that a in this particular case templated functions could be sufficient but to underline how to implement the most general case a class is used.

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-def-start
    :end-before: rep-options-def-end
    :dedent: 4

The specific implementation of the two options is done in with template specialization

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-impl-start
    :end-before: rep-options-impl-end
    :dedent: 4

The switch between the two behaviours is done in the :cpp:func:`compute() <rascal::CalculatorSortedCoulomb::compute()>` function which chooses the right type of computation method to compute the representation for each structure manager. At the moment the :cpp:func:`compute_loop() <rascal::CalculatorSortedCoulomb::compute_loop()>` function has to be copied to every new calculator to allow iterations over structure managers. It is not necessary to understand this code, if one is only interested to integegrate a new representation into rascal.

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: compute-loop-begin
    :end-before: compute-loop-end
    :dedent: 4


The actual implementation of the representation is in the function :cpp:func:`compute_impl() <rascal::CalculatorSortedCoulomb::compute_impl()>` witch the computation method as template argument.

.. literalinclude:: ../../../src/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-compute-start
    :end-before: rep-options-compute-end
    :dedent: 2

Write the python bindings of a new representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use ``pybind11`` to handle the generation of python bindings. To make a new representation available for the use on the python side, you need to explicitly register your representation calculator for every possible structure manager stack that you want to make available in the :cpp:func:`add_representation_calculators() <add_representation_calculators>` function in ``bindings/bind_py_representation_calculator.cc``. Here is an example on how it is done for the sorted coulomb representation.

.. literalinclude:: ../../../bindings/bind_py_representation_calculator.cc
    :language: c++
    :start-after: rep-bind-start
    :end-before: rep-bind-end
    :dedent: 2

The last step is to write a python class in ``bindings/rascal/representations/`` to simplify the use of the representation from the python side. You can use the implentation of the :class:`SortedCoulombMatrix` in the ``coulomb_matrix.py`` as a template. 



.. [#one] http://www.qmlcode.org/qml.html#module-qml.representations.

.. [#two] Rupp, M., Tkatchenko, A., Müller, K.-R., & von Lilienfeld, O. A. (2011).
        Fast and Accurate Modeling of Molecular Atomization Energies with Machine Learning.
        Physical Review Letters, 108(5), 58301. https://doi.org/10.1103/PhysRevLett.108.058301
