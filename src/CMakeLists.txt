


set(_proteus_sources
        basic_types.cc
        )
set(_proteus_headers
        basic_types.h
        )


pybind11_add_module(proteus module.cc  ${_proteus_sources} )
