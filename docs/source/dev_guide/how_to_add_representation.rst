.. _how_to_add_representation:

.. contents::
   :local:

How to add a new representation
-------------------------------



Write a RepresentationManager
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A representation manager is an object that builds a representation of the atomic structure contained in *a* structure manager. This class:

- inherits publicly from :cpp:class:`RepresentationManagerBase <rascal::RepresentationManagerBase>` to follow its interface and use some of the common utilies shared by such class.

- is templated by a structure manager to be able to build the representation efficiently.

- uses the :cpp:class:`Property <rascal::Property>` class to store the representation's features.

The behaviour of a representation manager is defined at construction by a structure manager and a structure containing a set of parameters. Within these parameters some
  - define behaviours or options for the class, e.g. choosing a particular implementation from a set of conceptually equivalent methods.
  - set the hyperparameters of the representation.

Note that there is one representation manager per structure manager.

To illustrate the basic structure that a new representation that would be implemented in ``representation_manager_custom.hh`` should follow, let's take the example of the :cpp:class:`RepresentationManagerSortedCoulomb <rascal::RepresentationManagerSortedCoulomb>`. A detailed discussion of the sorted coulomb matrix representation can be found in [#one]_ and [#two]_.


The representation starts with the definition of some useful short hands

.. literalinclude:: ../../../src/representations/representation_manager_sorted_coulomb.hh
    :language: c++
    :start-after: rep-preamble-start
    :end-before: rep-preamble-end
    :dedent: 2

followed by the definition of its constructors and destructor

.. literalinclude:: ../../../src/representations/representation_manager_sorted_coulomb.hh
    :language: c++
    :start-after: rep-construc-start
    :end-before: rep-construc-end
    :dedent: 2

the declaration of the concrete implementation of the RepresentationManager interface

.. literalinclude:: ../../../src/representations/representation_manager_sorted_coulomb.hh
    :language: c++
    :start-after: rep-interface-start
    :end-before: rep-interface-end
    :dedent: 2

and the declaration of some functions for internal use in the protected section.

The end of the class contains the different internal variables needed by the class

.. literalinclude:: ../../../src/representations/representation_manager_sorted_coulomb.hh
    :language: c++
    :start-after: rep-variables-start
    :end-before: rep-variables-end
    :dedent: 2


Write the python bindings of the new representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use ``pybind11`` to handle the generation of python bindings. The first step to register the new representation is to include the file in ``bindings/bind_include.hh``. Then, you need to explicitly register your representation manager for every possible structure manager stack that you want to make available to the python side in the :cpp:func:`add_representation_managers() <add_representation_managers>`. Here is an example on how it is done for the sorted coulomb representation

.. literalinclude:: ../../../bindings/bind_py_representation_manager.cc
    :language: c++
    :start-after: rep-bind-start
    :end-before: rep-bind-end
    :dedent: 2

The last step is to write a python class in ``bindings/rascal/representations/`` to simplify the use of the representation from the python side. You can use :class:`SortedCoulombMatrix` as a template



Efficient Implementation of Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The switch between several implementations of conceptually equivalent parts of a representation can be implemented through several mechanisms such a virtual inheritance. We detail here how to implement such switch efficiently using the :cpp:class:`RepresentationManagerSortedCoulomb <rascal::RepresentationManagerSortedCoulomb>` as an example.

To make the coulomb matrix invariant with respect to permutation Ref. [#two]_ proposes to sort the upper triangular part of the coulomb matrix according to the norm of each rows or the distance from the central atom (see [#one]_ for details).

The implementation of these two behaviour is encapsulated in the :cpp:class:`SortCoulomMatrix <rascal::internal::SortCoulomMatrix>` class and the choice between them is done with a template parameter using template specialization. Note that a in this particular case templated functions could be sufficient but to underline how to implement the most general case a class is used.

.. literalinclude:: ../../../src/representations/representation_manager_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-def-start
    :end-before: rep-options-def-end
    :dedent: 4

The specific implementation of the two options is done in with template specialization

.. literalinclude:: ../../../src/representations/representation_manager_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-impl-start
    :end-before: rep-options-impl-end
    :dedent: 4

Finally, the switch between the two behaviours is done in the :cpp:func:`compute() <rascal::RepresentationManagerSortedCoulomb::compute()>` by calling the templated function :cpp:func:`compute_helper() <rascal::RepresentationManagerSortedCoulomb::compute_helper()>` where the computation of the representation is actually implemented

.. literalinclude:: ../../../src/representations/representation_manager_sorted_coulomb.hh
    :language: c++
    :start-after: rep-options-compute-start
    :end-before: rep-options-compute-end
    :dedent: 2



.. [#one] http://www.qmlcode.org/qml.html#module-qml.representations.

.. [#two] Rupp, M., Tkatchenko, A., Müller, K.-R., & von Lilienfeld, O. A. (2011).
        Fast and Accurate Modeling of Molecular Atomization Energies with Machine Learning.
        Physical Review Letters, 108(5), 58301. https://doi.org/10.1103/PhysRevLett.108.058301
