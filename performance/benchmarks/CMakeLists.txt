file( GLOB binaries "${CMAKE_SOURCE_DIR}/performance/benchmarks/*.cc")
list( APPEND binaries "${CMAKE_SOURCE_DIR}/docs/source/tutorials/benchmark_tutorial.cc")
foreach(binaryfile ${binaries})
  get_filename_component(binaryname  ${binaryfile} NAME_WE)
  add_executable(${binaryname} ${binaryfile})
  target_link_libraries(${binaryname} "${LIBRASCAL_NAME}")
  target_link_libraries(${binaryname} Threads::Threads)
  target_link_libraries(${binaryname} benchmark)

  # surpresses the "'CSVReporter' is deprecated" errors when compiling benchmark.cc
  set_target_properties(${binaryname} PROPERTIES COMPILE_OPTIONS "-Wno-deprecated-declarations" )
endforeach(binaryfile ${binaries})
