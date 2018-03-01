set (PY_BINDING_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_module.cc
  )

pybind11_add_module(pyProteus ${PY_BINDING_SRCS})
target_link_libraries(pyProteus PRIVATE proteus)
# Want to rename the output, so that the python module is called _proteus
set_target_properties(pyProteus PROPERTIES OUTPUT_NAME _proteus)
