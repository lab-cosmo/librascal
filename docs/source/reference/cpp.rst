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

Kernels
~~~~~~~

Kernel
^^^^^^
 .. doxygenclass:: rascal::Kernel
    :project: rascal
    :members: 
    
 .. doxygenstruct:: rascal::internal::KernelImpl< internal::KernelType::Cosine >
    :project: rascal
    :members:
    

SparseKernel
^^^^^^^^^^^^

 .. doxygenclass:: rascal::SparseKernel
    :project: rascal
    :members: 


 .. doxygenstruct:: rascal::internal::SparseKernelImpl< internal::SparseKernelType::GAP >
    :project: rascal
    :members:

Math utilities (namespace :cpp:class:`rascal::math`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .. doxygenclass:: rascal::math::SphericalHarmonics
    :project: rascal
    :members:
    :private-members:

 .. doxygenclass:: rascal::math::Hyp1f1
    :project: rascal
    :members:

Index
~~~~~

 .. doxygennamespace:: rascal
    :project: rascal
    :members:
    :outline:
    :no-link:

