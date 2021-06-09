
file(GLOB examples "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

foreach(_file_ ${examples})
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_executable(${_name_} ${_file_})
    target_link_libraries(${_name_} "${LIBRASCAL_NAME}")
endforeach()

# BEGIN for libtorch
file(GLOB examples "${CMAKE_CURRENT_SOURCE_DIR}/torch/*.cc")
foreach(_file_ ${examples})
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_executable(${_name_} ${_file_})
    target_link_libraries(${_name_} "${LIBRASCAL_NAME}")
    target_link_libraries(${_name_} "$ENV{TORCH_LIBRARIES}")
    #set_property(TARGET ${_name_} PROPERTY CXX_STANDARD 14)
endforeach()
# END for libtorch

configure_file("${CMAKE_SOURCE_DIR}/reference_data/inputs/crystal_structure.json" ${CMAKE_BINARY_DIR} COPYONLY)
configure_file("${CMAKE_SOURCE_DIR}/reference_data/inputs/alanine-X.json" ${CMAKE_BINARY_DIR} COPYONLY)
configure_file("${CMAKE_SOURCE_DIR}/reference_data/inputs/alanine-center-select.json" ${CMAKE_BINARY_DIR} COPYONLY)
