function(_add_boost_lib)
  set(options )
  set(oneValueArgs NAME)
  set(multiValueArgs SOURCES LINK DEFINE DEFINE_PRIVATE CXXFLAGS_PRIVATE INCLUDE_PRIVATE)
  cmake_parse_arguments(BOOSTLIB "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN})
  add_library(Boost_${BOOSTLIB_NAME} STATIC ${BOOSTLIB_SOURCES})
  add_library(Boost::${BOOSTLIB_NAME} ALIAS Boost_${BOOSTLIB_NAME})
  set_target_properties(Boost_${BOOSTLIB_NAME} PROPERTIES
    OUTPUT_NAME "boost_${BOOSTLIB_NAME}"
    FOLDER "Boost"
  )
  if(NOT BOOST_STANDALONE)
    set_target_properties(Boost_${BOOSTLIB_NAME} PROPERTIES EXCLUDE_FROM_ALL 1)
  endif()
  target_link_libraries(Boost_${BOOSTLIB_NAME} PUBLIC Boost::unit_test_framework)
  if(MSVC)
    target_compile_options(Boost_${BOOSTLIB_NAME} PRIVATE /w)
  else()
    target_compile_options(Boost_${BOOSTLIB_NAME} PRIVATE -w)
  endif()
  if(BOOSTLIB_LINK)
    target_link_libraries(Boost_${BOOSTLIB_NAME} PUBLIC ${BOOSTLIB_LINK})
  endif()
  if(BOOSTLIB_DEFINE)
    target_compile_definitions(Boost_${BOOSTLIB_NAME} PUBLIC ${BOOSTLIB_DEFINE})
  endif()
  if(BOOSTLIB_DEFINE_PRIVATE)
    target_compile_definitions(Boost_${BOOSTLIB_NAME} PRIVATE ${BOOSTLIB_DEFINE_PRIVATE})
  endif()
  if(BOOSTLIB_CXXFLAGS_PRIVATE)
    target_compile_options(Boost_${BOOSTLIB_NAME} PRIVATE ${BOOSTLIB_CXXFLAGS_PRIVATE})
  endif()
  if(BOOSTLIB_INCLUDE_PRIVATE)
    target_include_directories(Boost_${BOOSTLIB_NAME} PRIVATE ${BOOSTLIB_INCLUDE_PRIVATE})
  endif()
endfunction()


_add_boost_lib(
  NAME unit_test_framework
  SOURCES
    ${_boost_test_external_dir}/src/compiler_log_formatter.cpp
    ${_boost_test_external_dir}/src/debug.cpp
    ${_boost_test_external_dir}/src/decorator.cpp
    ${_boost_test_external_dir}/src/execution_monitor.cpp
    ${_boost_test_external_dir}/src/framework.cpp
    ${_boost_test_external_dir}/src/plain_report_formatter.cpp
    ${_boost_test_external_dir}/src/progress_monitor.cpp
    ${_boost_test_external_dir}/src/results_collector.cpp
    ${_boost_test_external_dir}/src/results_reporter.cpp
    ${_boost_test_external_dir}/src/test_tools.cpp
    ${_boost_test_external_dir}/src/test_tree.cpp
    ${_boost_test_external_dir}/src/unit_test_log.cpp
    ${_boost_test_external_dir}/src/unit_test_main.cpp
    ${_boost_test_external_dir}/src/unit_test_monitor.cpp
    ${_boost_test_external_dir}/src/unit_test_parameters.cpp
    ${_boost_test_external_dir}/src/junit_log_formatter.cpp
    ${_boost_test_external_dir}/src/xml_log_formatter.cpp
    ${_boost_test_external_dir}/src/xml_report_formatter.cpp
)