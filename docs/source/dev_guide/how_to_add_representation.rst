.. _how_to_add_representation:

.. contents::
   :local:

How to add a new representation
-------------------------------



Write a Calculator for a representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A calculator of a representation is an object that builds a representation of the atomic structure contained in *a* structure manager. This class:

- inherits publicly from :cpp:class:`CalculatorBase <rascal::CalculatorBase>` to follow its interface and use some of the common utilies shared by such class.

- uses a :cpp:func:`compute() <compute()>` function expecting one or a list of structure manager(s) to build the representation efficiently and attach the resulting features to their respective structure manager.

The behaviour of a representation calculator is defined at construction by hyperparameters which are specific to the representation.
These hyperparameters can:

- define a particular implementation from a set of conceptually equivalent methods, implemented using primitive enums

- be scalar values or on/off features used as hyperparameters in the calculation of the represenation, implemented using primitive types


To illustrate the basic structure that a new representation that would be implemented in ``calculator_representation_name.hh`` should follow, let's take the example of the :cpp:class:`CalculatorSortedCoulomb <rascal::CalculatorSortedCoulomb>`. A detailed discussion of the sorted coulomb matrix representation can be found in Ref. [#one]_ and Ref. [#two]_.


The representation starts with the definition of some useful short hands

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-preamble-start
    :end-before: rep-preamble-end
    :dedent: 2

followed by the definition of its constructors and destructor

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-construc-start
    :end-before: rep-construc-end
    :dedent: 4

the declaration of the concrete implementation of the calculator interface

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-interface-start
    :end-before: rep-interface-end
    :dedent: 4

and the declaration of some functions for internal use in the protected section. The end of the class contains the different internal variables needed by the class

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-variables-start
    :end-before: rep-variables-end
    :dedent: 4

Implementing several types of the same representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Often there are different ways to implement the same fundamental representation. To be able to do this in rascal, the compute function is used with the type of the implementation as template argument (or multiple if representation requires this). The switch between several implementations of conceptually equivalent parts of a representation can be implemented through several mechanisms such a virtual inheritance. We detail here how to implement such switch efficiently using the :cpp:class:`CalculatorSortedCoulomb <rascal::CalculatorSortedCoulomb>` as an example.

To make the coulomb matrix invariant with respect to permutation Ref. [#two]_ proposes to sort the upper triangular part of the coulomb matrix according to the norm of each rows or the distance from the central atom (see [#one]_ for details).

The implementation of these two behaviour is encapsulated in the :cpp:class:`SortCoulomMatrix <rascal::internal::SortCoulomMatrix>` class and the choice between them is done with a template parameter using template specialization. Note that in this particular case templated functions could be sufficient but to underline how to implement the most general case a class is used.

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-def-start
    :end-before: rep-options-def-end
    :dedent: 4

The specific implementation of the two options is done in with template specialization

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-impl-start
    :end-before: rep-options-impl-end
    :dedent: 4

The switch between the two behaviours is done in the :cpp:func:`compute() <rascal::CalculatorSortedCoulomb::compute()>` function which chooses the right type of computation method to compute the representation for each structure manager. At the moment the :cpp:func:`compute_loop() <rascal::CalculatorSortedCoulomb::compute_loop()>` function has to be copied to every new calculator to allow iterations over structure managers. It is not necessary to understand this code, if one is only interested to integegrate a new representation into rascal.

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: compute-loop-begin
    :end-before: compute-loop-end
    :dedent: 4


The actual implementation of the representation is in the function :cpp:func:`compute_impl() <rascal::CalculatorSortedCoulomb::compute_impl()>` which is templated with the computation method.

.. literalinclude:: ../../../src/rascal/representations/calculator_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-compute-start
    :end-before: rep-options-compute-end
    :dedent: 2

Write the python bindings of a new representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use ``pybind11`` to handle the generation of the python bindings. To make a new representation available to the python users, you need to explicitly register your representation calculator. Since the function :cpp:func:`compute() <rascal::CalculatorSortedCoulomb::compute()>` is templated with the type of the input structure manager, it has to be binded explicitly for every possible structure manager stack type that you want to allow as input. To do this you need to add you own binding code in the :cpp:func:`add_representation_calculators() <add_representation_calculators>` function in ``bindings/bind_py_representation_calculator.cc`` using the available infrastructure. Here is an example on how it is done for the sorted coulomb representation.

.. literalinclude:: ../../../bindings/bind_py_representation_calculator.cc
    :language: c++
    :start-after: rep-bind-start
    :end-before: rep-bind-end
    :dedent: 4


.. Note::

    To be able to use a particular structure manager stack in python, it also has to be binded. In the case a valid structure manager stack for your representation is not already binded, you will have to register it too in the :cpp:func:`add_structure_managers() <add_structure_managers>` function in ``bindings/bind_py_structure_manager.cc`` like in this example:
    
    .. literalinclude:: ../../../bindings/bind_py_structure_manager.cc
        :language: c++
        :start-after: struc-man-bind-start
        :end-before: struc-man-bind-end
        :dedent: 4


The last step is to write a python class in ``bindings/rascal/representations/`` to simplify the use of the representation from the python side. You can use the implentation of the :class:`SortedCoulombMatrix` in ``coulomb_matrix.py`` as a template.


Implement and test gradients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In principle, all libRascal representations should implement gradients with
respect to the atomic positions.  Currently the only representations to do so
are the :cpp:class:`SphericalExpansion` and :cpp:class:`SphericalInvariants`
(tensor order 0, body order 0--1 a.k.a. "RadialSpectrum" and "PowerSpectrum") in
the GTO and DVR radial basis.  Until we come up with a general, standard way of
implementing gradients for any representation, please see those implementations
for guidance.

Once you've implemented the gradient -- or derivative -- of any function in
libRascal, you must test that it actually corresponds to the gradient of the
function that it purpots to be.  A finite-difference testing function is
provided for this purpose; see :ref:`testing_gradients` for details.


.. [#one] http://www.qmlcode.org/qml.html#module-qml.representations.

.. [#two] Rupp, M., Tkatchenko, A., Müller, K.-R., & von Lilienfeld, O. A. (2011).
        Fast and Accurate Modeling of Molecular Atomization Energies with Machine Learning.
        Physical Review Letters, 108(5), 58301. https://doi.org/10.1103/PhysRevLett.108.058301
