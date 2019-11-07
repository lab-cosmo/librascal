find_package(Boost REQUIRED COMPONENTS unit_test_framework )

option(RASCAL_TESTS_USE_VALGRIND "run all C++ tests within valgrind" OFF)

if (${RASCAL_TESTS_USE_VALGRIND})
    find_program(VALGRIND_EXECUTABLE "valgrind")
    if ("${VALGRIND_EXECUTABLE}" STREQUAL "VALGRIND_EXECUTABLE-NOTFOUND")
        message(FATAL_ERROR "Could not find valgrind")
    endif()
    set(TEST_RUNNER "${VALGRIND_EXECUTABLE}"
        # fail the test if valgrind fails
        "--error-exitcode=100"
        # check for all kinds of leaks
        "--leak-check=full"
        "--show-leak-kinds=all"
        "--errors-for-leak-kinds=all"
        # track leaks origin to improve error messages
        "--track-origins=yes"
        # use dsymutil to extract debug symbols on macOS
        "--dsymutil=yes"
        # remove errors from the OS
        "--suppressions=${CMAKE_CURRENT_SOURCE_DIR}/valgrind.suppressions"
    )

else()
    set(TEST_RUNNER "")
endif()

add_library(rascal_tests_main STATIC helpers/tests_main.cc)
target_compile_definitions(rascal_tests_main PUBLIC -DBOOST_TEST_DYN_LINK)
target_link_libraries(rascal_tests_main ${Boost_LIBRARIES} "${LIBRASCAL_NAME}")

file(GLOB tests "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cc")

foreach(_file_ ${tests})
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_executable(${_name_} ${_file_})
    target_link_libraries(${_name_} rascal_tests_main)
    target_compile_options(${_name_} PRIVATE -Werror)

    add_test(
        NAME ${_name_}
        COMMAND ${TEST_RUNNER} $<TARGET_FILE:${_name_}> -- --report_level=detailed --build_info=TRUE
        # Run the tests in the root directory to give them access to the
        # reference data
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
endforeach()

if(BUILD_BINDINGS)
    add_test(
        NAME python_binding_test
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/python_binding_tests.py
        # Run the tests in the root directory to give them access to the
        # reference data
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    # Tell python where to find the rascal module
    set_tests_properties(python_binding_test
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}"
    )
endif()
