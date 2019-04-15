.. _property:

Property
--------

An object of the ``Property`` class is a container for a physical property assigned to a cluster. Only physical properties of the form of a scalar, vector or matrix can be expressed. As the property is a member of an adaptor, it is together with the adaptor bound to the same layer. This will be explained in greater detail in :ref:`how to add a new adaptor<add-property-to-adaptor>`. An example of how to declare different kinds of properties:

.. literalinclude:: ../../../examples/example_property.cc
    :language: c++
    :start-after: property-typedef-start
    :end-before: property-typedef-end
    :caption: examples/example_property.cc

As written in :ref:`the implementation of structure managers <structure-manager>`, an adaptor with the functionality to calculate a certain physical property has a ``Property`` object describing this property and making the property values available for access. A property can be constructed by handing the corresponding structure manager to the property and can then be filled with values by using the `push_back` member function as seen in the code below.

.. literalinclude:: ../../../examples/example_property.cc
    :language: c++
    :start-after: property-construction-start
    :end-before: property-construction-end
    :caption: examples/example_property.cc
    :dedent: 2

The values of the property can be accessed with a ``AtomRef`` or ``ClusterRef`` object.

.. literalinclude:: ../../../examples/example_property.cc
    :language: c++
    :start-after: property-access-start
    :end-before: property-access-end
    :caption: examples/example_property.cc
    :dedent: 2

Dynamic properties
******************

For properties with a dynamic size per cluster, the value ``Eigen::Dynamic`` should be used in `NbCol` and `NbRow`. When using a dynamic sized property one should be aware that dynamic properties are less performant than static ones.

TypedProperty
*************

For properties with unknown size (unknown `NbCol` and `NbRow`) at compilation time or any other use case where the size cannot be handed to the type definiton of the corresponding property, one can use ``TypedProperty``, where the assignment of the size happens in the constructor. 


.. literalinclude:: ../../../src/structure_managers/property_typed.hh
    :language: c++
    :start-after: typed-property-constructor-start
    :end-before: typed-property-constructor-end
    :caption: src/structure_managers/property_typed.hh
    :dedent: 4
