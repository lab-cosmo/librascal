download_external_project(BoostTest
  URL "https://github.com/boostorg/test.git"
  TAG "boost-${_BoostTest_version}"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${_BoostTest_external_dir}
  NO_UPDATE
  )


set(BOOST_TEST_ROOT_DIR ${_BoostTest_external_dir}/BoostTest)
get_filename_component(BOOST_TEST_ROOT_DIR_ABS ${_BoostTest_external_dir}/BoostTest ABSOLUTE)

# include globs
file(GLOB_RECURSE
     BOOST_UTF_HEADERS
     ${BOOST_TEST_ROOT_DIR}/include/*.hpp
     ${BOOST_TEST_ROOT_DIR}/include/*.ipp)

# organize files
foreach(_h IN LISTS BOOST_UTF_HEADERS)
  get_filename_component(_hh ${_h} ABSOLUTE)
  file(RELATIVE_PATH _v ${BOOST_TEST_ROOT_DIR_ABS}/include/boost/test ${_hh})
  get_filename_component(_v "${_v}" DIRECTORY)
  string(REPLACE "/" "\\" _v "${_v}")
  source_group(${_v} FILES ${_h})
endforeach()

set(BOOST_UTF_SRC
    # list specific to 1.58.0
    ${BOOST_TEST_ROOT_DIR}/src/compiler_log_formatter.cpp
    ${BOOST_TEST_ROOT_DIR}/src/debug.cpp
    ${BOOST_TEST_ROOT_DIR}/src/exception_safety.cpp
    ${BOOST_TEST_ROOT_DIR}/src/execution_monitor.cpp
    ${BOOST_TEST_ROOT_DIR}/src/framework.cpp
    ${BOOST_TEST_ROOT_DIR}/src/interaction_based.cpp
    ${BOOST_TEST_ROOT_DIR}/src/logged_expectations.cpp
    ${BOOST_TEST_ROOT_DIR}/src/plain_report_formatter.cpp
    ${BOOST_TEST_ROOT_DIR}/src/progress_monitor.cpp
    ${BOOST_TEST_ROOT_DIR}/src/results_collector.cpp
    ${BOOST_TEST_ROOT_DIR}/src/results_reporter.cpp
    ${BOOST_TEST_ROOT_DIR}/src/test_tools.cpp
    ${BOOST_TEST_ROOT_DIR}/src/unit_test_log.cpp
    ${BOOST_TEST_ROOT_DIR}/src/unit_test_main.cpp
    ${BOOST_TEST_ROOT_DIR}/src/unit_test_monitor.cpp
    ${BOOST_TEST_ROOT_DIR}/src/unit_test_parameters.cpp
    ${BOOST_TEST_ROOT_DIR}/src/unit_test_suite.cpp
    ${BOOST_TEST_ROOT_DIR}/src/xml_log_formatter.cpp
    ${BOOST_TEST_ROOT_DIR}/src/xml_report_formatter.cpp
    )

add_library(boost_test_framework STATIC ${BOOST_UTF_HEADERS} ${BOOST_UTF_SRC})

target_compile_definitions(boost_test_framework PUBLIC "-DBOOST_TEST_DYN_LINK=0")

if (SUPPRESS_LIBRARY_WARNINGS)
  target_include_directories(boost_test_framework SYSTEM PUBLIC ${BOOST_TEST_ROOT_DIR}/include/)
  message(STATUS "Suppressing warnings from library Boost Unit Test Framework")
else (SUPPRESS_LIBRARY_WARNINGS)
  target_include_directories(boost_test_framework PUBLIC ${BOOST_TEST_ROOT_DIR}/include/)
endif (SUPPRESS_LIBRARY_WARNINGS)

set(BOOST_TEST_INCLUDE_DIR ${BOOST_TEST_ROOT_DIR}/include/)

set_target_properties(boost_test_framework PROPERTIES FOLDER "UTF")

message(STATUS "Boost Unit Test Framework " ${_BoostTest_version})

