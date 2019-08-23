#message(WARNING "processed: ${_benchmark_external_dir}")
#TODO(alex) reinsert download_external_project
# If project does not exist 
#download_external_project(benchmark
#  URL "https://github.com/google/benchmark.git"
#  BACKEND GIT
#  THIRD_PARTY_SRC_DIR ${_benchmark_external_dir}
#  BUILD_COMMAND ""
#  ${_benchmark_UPDATE}
#  )

set(benchmark_SOURCE_DIR "${_benchmark_external_dir}/${package}")
#message(WARNING "processed: ${benchmark_SOURCE_DIR}")

# We do not require to run the tests of benchmark, but in case
#download_external_project(googletest
#  URL "https://github.com/google/googletest.git"
#  BACKEND GIT
#  THIRD_PARTY_SRC_DIR ${benchmark_SOURCE_DIR}
#  ${_googletest_update}
#  )

# If enabled the -std=c++03 flag is enabled when compiling cxx03_test.cc thus resulting in an error. We do not need the test anyway, but the cmake file does things I do not get outside this project, it is hard to locate the reason for any compiling error.
# TODO(all) enable testing is still an compitle option, but will not work because of the abovementioned error. Should I fix the compile error or rather make disable option. I am not sure how to do neither of them.
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable benchmark self-testing?")
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Enable benchmark install option?")
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "")


# BUILD_COMMAND "mkdir -R ${benchmark_SOURCE_DIR}/build && cmake -DCMAKE_BUILD_TYPE=Release --build ${benchmark_SOURCE_DIR}/build && make -C ${benchmark_SOURCE_DIR}/build"


add_subdirectory("${benchmark_SOURCE_DIR}")
include_directories("${benchmark_SOURCE_DIR}/include")

# surpresses the "'CSVReporter' is deprecated" errors when compiling benchmark.cc
set_target_properties(
    benchmark PROPERTIES
    COMPILE_FLAGS
      "-Wno-deprecated-declarations"
    )

#####################################

#include(ConfigSafeGuards)
#set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON CACHE BOOL "Download dependencies?")
#if (CMAKE_VERSION VERSION_LESS 3.2)
#    set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
#else()
#    set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
#endif()
#
#include(DownloadProject)
#download_project(PROJ                googlebenchmark
#                 GIT_REPOSITORY      https://github.com/google/benchmark.git
#                 GIT_TAG             v1.5.0
#                 ${UPDATE_DISCONNECTED_IF_AVAILABLE}
#)
#
#add_subdirectory(${googlebenchmark_SOURCE_DIR} ${googlebenchmark_BINARY_DIR})
#
#include_directories("${googlebenchmark_SOURCE_DIR}/include")
#
#include_directories(SYSTEM BEFORE
##${benchmark_SOURCE_DIR}
#${benchmark_SOURCE_DIR}/include)
##${benchmark_SOURCE_DIR}/build/src)

####################################

#if((CMAKE_CXX_COMPILER_ID MATCHES GNU) OR (CMAKE_CXX_COMPILER_ID MATCHES Clang))
#   set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11")
#   set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
#   set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
#endif()

#message(WARNING ${CMAKE_CXX_FLAGS})
#message(WARNING "${CMAKE_SOURCE_DIR}/build/build/external/benchmark/")
# I think this is done by our download_external_project function
#add_subdirectory(${benchmark_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/build/build/external/benchmark/)

#add_library(benchmark STATIC IMPORTED)
#set_property(TARGET benchmark PROPERTY IMPORTED_LOCATION ${benchmark_SOURCE_DIR}/build/src/libbenchmark.a)
#set_property(TARGET benchmark PROPERTY IMPORTED_LOCATION ${benchmark_SOURCE_DIR}/build/src/libbenchmark_main.a)
#add_dependencies("${WIGXJPF_NAME}" wigxjpf)

#mark_as_advanced_prefix(benchmark)
#mark_as_advanced(USE_PYTHON_INCLUDE_DIR)

set(${package}_FOUND TRUE CACHE INTERNAL "To avoid cyclic search" FORCE)
set(${package}_FOUND_EXTERNAL TRUE CACHE INTERNAL "" FORCE)
