.. _auto_cpp:

C++ documentation
-----------------

This list is incomplete. You can help by *expanding it*!

Representations
~~~~~~~~~~~~~~~

Spherical Expansion
^^^^^^^^^^^^^^^^^^^

 .. doxygenstruct:: rascal::internal::RadialContribution< RadialBasisType::GTO >
    :project: rascal
    :members:

 .. doxygenstruct:: rascal::internal::RadialContribution< RadialBasisType::DVR >
    :project: rascal
    :members:

 .. doxygenclass:: rascal::CalculatorSphericalExpansion
    :project: rascal
    :members:

Spherical Invariants
^^^^^^^^^^^^^^^^^^^^

 .. doxygenclass:: rascal::CalculatorSphericalInvariants
    :project: rascal
    :members:

Spherical Covariants
^^^^^^^^^^^^^^^^^^^^

 .. doxygenclass:: rascal::CalculatorSphericalCovariants
    :project: rascal
    :members:

Math utilities (namespace :cpp:class:`rascal::math`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .. doxygenclass:: rascal::math::SphericalHarmonics
    :project: rascal
    :members:

 .. doxygenclass:: rascal::math::Hyp1f1
    :project: rascal
    :members:

Gradient tests (developer)
^^^^^^^^^^^^^^^^^^^^^^^^^^

All representations in libRascal should implement, or plan to implement,
gradients (derivatives w.r.t. the Cartesian positions of all the atoms) so that
they can be used to run dynamics.  This is often a complex and error-prone task,
so a finite-difference gradient checker is provided to check the gradients of
any representation calculator -- or any math function in general -- and ensure
that the analytical and finite-difference gradients match up.

To check the gradient of a new representation calculator, it should suffice to
use the classes RepresentationManagerGradientCalculator (to provide the function
and its gradient) and RepresentationManagerGradientFixture (to assist in
iterating over the atoms of the structure).  An example of its usage is shown
below, excerpted from :file:`tests/test_calculator.cc`:

.. literalinclude:: ../../../tests/test_calculator.cc
   :language: cpp
   :lines: 430-436, 441-444
   :dedent: 2


Index
~~~~~

 .. doxygennamespace:: rascal
    :project: rascal
    :members:

