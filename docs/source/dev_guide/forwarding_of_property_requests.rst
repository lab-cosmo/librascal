.. _forwarding_of_property_request:

.. contents::
   :local:

Forwarding of property requests
-------------------------------

This c++ file is a simplification of our CRTP structure to show how the forwarding of requests through the manager stack works. The forwarding mechanism is used when getting a property with the `get_property` function and for the `get_distance` and `get_direction_vector` function.

.. literalinclude:: ../../../examples/forwarding_of_property_requests.cc
    :language: c++
