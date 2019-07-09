download_external_project(googlebenchmark
  URL "https://github.com/google/benchmark.git"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${_googlebenchmark_external_dir}
  ${_googlebenchmark_update}
  )

set(SOURCE_DIR "${_googlebenchmark_external_dir}/${package}")

download_external_project(googletest
  URL "https://github.com/google/googletest.git"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${SOURCE_DIR}
  BUILD_COMMAND "mkdir -R ${SOURCE_DIR}/build && cmake -DCMAKE_BUILD_TYPE=Release --build ${SOURCE_DIR}/build && make -C ${SOURCE_DIR}/build"
  ${_googletest_update}
  )

#add_subdirectory(${_googlebenchmark_external_dir}/googlebenchmark)

include_directories(SYSTEM ${benchmark_INCLUDE_DIRS})
#file(READ "${_googlebenchmark_external_dir}/googlebenchmark/Eigen/src/Core/util/Macros.h" _eigen_version_header)

#add_library("${GOOGLEBENCHMARK_NAME}" STATIC IMPORTED)
#set_property(TARGET "${GOOGLEBENCHMARK_NAME}" PROPERTY IMPORTED_LOCATION ${SOURCE_DIR}/build/src/libbenchmark.a)
#set_property(TARGET "${GOOGLEBENCHMARK_NAME}" PROPERTY IMPORTED_LOCATION ${SOURCE_DIR}/build/src/libbenchmark_main.a)
#add_dependencies("${GOOGLEBENCHMARK_NAME}" ${package})


set(GOOGLEBENCHMARK_INCLUDE_DIR ${_googlebenchmark_external_dir}/googlebenchmark/include CACHE PATH "Googlebenchmark include directory")
#target_include_directories(googlebenchmark SYSTEM INTERFACE
#  $<BUILD_INTERFACE:${GOOGLEBENCHMARK_INCLUDE_DIR}>)

set(${package}_FOUND TRUE CACHE INTERNAL "To avoid cyclic search" FORCE)
set(${package}_FOUND_EXTERNAL TRUE CACHE INTERNAL "" FORCE)
