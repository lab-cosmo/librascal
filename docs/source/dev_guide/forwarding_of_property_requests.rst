.. _forwarding_of_property_request:

.. contents::
   :local:

Forwarding of property requests
-------------------------------

This c++ file is a simplification of our CRTP structure to show how the forwarding of requests through the manager stack work. The forwarding mechanism is used when getting a property with the `get_property_ptr` or `get_property_ref` function and for the `get_distance` and `get_direction_vector` function.

.. literalinclude:: forwarding_of_property_requests.cc
    :language: c++
