
file(GLOB examples "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

foreach(_file_ ${examples})
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_executable(${_name_} ${_file_})
    target_link_libraries(${_name_} "${LIBRASCAL_NAME}")
endforeach()

configure_file("${CMAKE_SOURCE_DIR}/examples/reference_data/crystal_structure.json" ${CMAKE_BINARY_DIR} COPYONLY)
configure_file("${CMAKE_SOURCE_DIR}/examples/reference_data/alanine-X.json" ${CMAKE_BINARY_DIR} COPYONLY)
configure_file("${CMAKE_SOURCE_DIR}/tests/reference_data/alanine-center-select.json" ${CMAKE_BINARY_DIR} COPYONLY)
