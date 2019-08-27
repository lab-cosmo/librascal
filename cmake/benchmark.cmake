# TODO(all) download external project always downloads this one, but other packages like Eigen are not always downloaded, does anyone have a clue?
download_external_project(benchmark
  URL "https://github.com/google/benchmark.git"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${_benchmark_external_dir}
  BUILD_COMMAND ""
  ${_benchmark_UPDATE}
  )

set(benchmark_SOURCE_DIR "${_benchmark_external_dir}/${package}")


# TODO(alex) fix BENCHMARK_ENABLE_TESTING problem:

# If we want to include tests then:
# If enabled the -std=c++03 flag is enabled when compiling cxx03_test.cc thus resulting in an error. We do not need the test anyway, but the cmake file does things I do not get outside this project, it is hard to locate the reason for any compiling error.

# We do not require to run the tests of benchmark, but in case
#download_external_project(googletest
#  URL "https://github.com/google/googletest.git"
#  BACKEND GIT
#  THIRD_PARTY_SRC_DIR ${benchmark_SOURCE_DIR}
#  ${_googletest_update}
#  )


# If we want to not include it, how to disable an option in cmake?

#set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable benchmark self-testing?")
#set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Enable benchmark install option?")
#set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "")



add_subdirectory("${benchmark_SOURCE_DIR}")
include_directories("${benchmark_SOURCE_DIR}/include")

# surpresses the "'CSVReporter' is deprecated" errors when compiling benchmark.cc
set_target_properties(
    benchmark PROPERTIES
    COMPILE_FLAGS
      "-Wno-deprecated-declarations"
    )

set(${package}_FOUND TRUE CACHE INTERNAL "To avoid cyclic search" FORCE)
set(${package}_FOUND_EXTERNAL TRUE CACHE INTERNAL "" FORCE)
