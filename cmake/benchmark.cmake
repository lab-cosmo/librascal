download_external_project(benchmark
    URL "https://github.com/google/benchmark.git"
    TAG "v${_benchmark_version}"
    BACKEND GIT
    THIRD_PARTY_SRC_DIR ${_benchmark_external_dir}
    BUILD_COMMAND ""
    ${_benchmark_UPDATE}
    )

set(benchmark_SOURCE_DIR "${_benchmark_external_dir}/${package}")

# turns BENCHMARK settings off which are not targeted
set(BENCHMARK_ENABLE_TESTING CACHE BOOL OFF)
set(BENCHMARK_ENABLE_GTEST_TESTS CACHE BOOL OFF)
set(BENCHMARK_ENABLE_INSTALL CACHE BOOL OFF)

mark_as_advanced(BENCHMARK_ENABLE_TESTING)
mark_as_advanced(BENCHMARK_BUILD_32_BITS)
mark_as_advanced(BENCHMARK_DOWNLOAD_DEPENDENCIES)
mark_as_advanced(BENCHMARK_ENABLE_ASSEMBLY_TESTS)
mark_as_advanced(BENCHMARK_ENABLE_EXCEPTIONS)
mark_as_advanced(BENCHMARK_ENABLE_GTEST_TESTS)
mark_as_advanced(BENCHMARK_ENABLE_INSTALL)
mark_as_advanced(BENCHMARK_ENABLE_LTO)
mark_as_advanced(BENCHMARK_USE_LIBCXX)

# HACK around https://github.com/google/benchmark/issues/773
# this might fail if we use another compiler than GCC/clang
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_subdirectory("${benchmark_SOURCE_DIR}")
include_directories(SYSTEM "${benchmark_SOURCE_DIR}/include")

# surpresses the "'CSVReporter' is deprecated" some warnings in google benchmark code
target_compile_options(benchmark PUBLIC -Wno-deprecated-declarations -Wno-effc++)

set(${package}_FOUND TRUE CACHE INTERNAL "To avoid cyclic search" FORCE)
set(${package}_FOUND_EXTERNAL TRUE CACHE INTERNAL "" FORCE)
