set(RUN_BENCHMARKS_FLAGS "" CACHE
    "Additional flags to run all benchmarks. See ./benchmark_interpolator"
    "--help for available flags")
mark_as_advanced(RUN_BENCHMARKS_FLAGS)

file(GLOB benchmarks "${CMAKE_SOURCE_DIR}/performance/benchmarks/*.cc")
list(APPEND benchmarks "${CMAKE_SOURCE_DIR}/docs/source/tutorials/benchmark_tutorial.cc")

set(ALL_CXX_BENCHMARKS "")
foreach(benchmark ${benchmarks})
  get_filename_component(name  ${benchmark} NAME_WE)
  add_executable(${name} ${benchmark})
  target_link_libraries(${name} "${LIBRASCAL_NAME}")
  target_link_libraries(${name} Threads::Threads)
  target_link_libraries(${name} benchmark)
  list(APPEND ALL_CXX_BENCHMARKS ${name})
endforeach()

set(CXX_BENCH_FLAGS "--benchmark_counters_tabular=true" ${RUN_BENCHMARKS_FLAGS})

add_custom_target(cpp_benchmarks
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/benchmark_interpolator ${CXX_BENCH_FLAGS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS ${ALL_CXX_BENCHMARKS}
)

add_custom_target(benchmarks DEPENDS cpp_benchmarks)

if (BUILD_BINDINGS)
    add_custom_target(python_benchmarks
        # Assume bash-like shell to pass PYTHONPATH
        COMMAND "PYTHONPATH=${CMAKE_BINARY_DIR}" ${PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/benchmark_spherical_representations.py
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        DEPENDS copy_python_bindings
    )

    add_dependencies(benchmarks python_benchmarks)
endif()
