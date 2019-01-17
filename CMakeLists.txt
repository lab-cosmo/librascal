# =============================================================================
# file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch
# @author Felix Musil <felix.musil@epfl.ch>
# @author Markus Stricker <markus.stricker@epfl.ch>
#
# @date   17 Jan 2019
#
# @brief  Main configuration file
#
# @section LICENSE
#
# Copyright © 2018 Till Junge
#
# Rascal is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# Rascal is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Rascal; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
# =============================================================================

cmake_minimum_required(VERSION 3.0)

project(rascal)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS ON)
set(RASCAL_PYTHON_MAJOR_VERSION 3)


add_compile_options(-Wall -Wextra)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
# to avoid conflicts with built-in functions
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-builtin")


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(rascalTools)


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
    OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  add_compile_options(-Wno-missing-braces)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  # add_compile_options(-Wno-non-virtual-dtor)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-non-virtual-dtor")
  string( TOLOWER "${CMAKE_BUILD_TYPE}" build_type )
  if ("release" STREQUAL "${build_type}" )
    add_compile_options(-march=native)
  endif()
  if ("debug" STREQUAL "${build_type}" )
    add_compile_options(-O0)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

add_external_package(Eigen3 VERSION 3.3.4 CONFIG)
add_external_package(pybind11 VERSION 2.2.2 CONFIG)

find_package(PythonInterp ${RASCAL_PYTHON_MAJOR_VERSION} REQUIRED)

# Try to find optional cpp linter
include(cpplint)

# Check for clang-format
include(clangformat)

include_directories(
  ${CMAKE_SOURCE_DIR}/tests
  ${CMAKE_SOURCE_DIR}/src
 )

include_directories(SYSTEM
  ${CMAKE_SOURCE_DIR}/external
  )

set(ENABLE_DOC "OFF" CACHE BOOL
  "If on, a target 'dev_doc' (which builds the documentation) is added")
set(BUILD_TESTS "ON" CACHE BOOL
  "If on, builds the unit tests")
set(SUPPRESS_LIBRARY_WARNINGS "ON" CACHE BOOL
  "Suppress warnings from  Boost library")

set(INSTALL "build" CACHE STRING
  "Use the default install path")

set(PYTHON_BUILD_PATH "${CMAKE_BINARY_DIR}/rascal/" CACHE STRING
    "Default install folder")


if (${BUILD_TESTS})

  find_package(Boost REQUIRED COMPONENTS unit_test_framework )
  enable_testing()

  ####################################################################
  # build c++ tests and add to test suite
  file( GLOB TEST_SRCS "${CMAKE_SOURCE_DIR}/tests/test_*.cc")

  add_executable(main_test_suite tests/main_test_suite.cc ${TEST_SRCS})

  # https://stackoverflow.com/questions/7781374/boost-test-linking
  ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

  target_include_directories(main_test_suite SYSTEM PRIVATE ${Boost_INCLUDE_DIR})
  target_link_libraries(main_test_suite ${Boost_LIBRARIES} rascal)
  add_test(main_test_suite main_test_suite
    --report_level=detailed
    --build_info=TRUE
    )

  # copy reference data
  configure_file("tests/alanine-X.json" ${CMAKE_BINARY_DIR} COPYONLY)
  configure_file("tests/simple_cubic_9.json" ${CMAKE_BINARY_DIR} COPYONLY)
  configure_file("examples/crystal_structure.json" ${CMAKE_BINARY_DIR} COPYONLY)
  file(COPY "tests/reference_data" DESTINATION ${CMAKE_BINARY_DIR})

  ####################################################################
  # copy and add python test to suite
  file( GLOB PY_TEST_SRCS "${CMAKE_SOURCE_DIR}/tests/python_*.py")

  foreach(pytest ${PY_TEST_SRCS})
    get_filename_component(pytest_name ${pytest} NAME)
    configure_file(
      ${pytest}
      "${CMAKE_BINARY_DIR}/${pytest_name}"
      COPYONLY)
  endforeach(pytest ${PY_TEST_SRCS})
  add_test(python_binding_test python_binding_tests.py)

endif(${BUILD_TESTS})

# compile the library
# add_compile_options(-Werror) # warnings are errors!
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/src/
  )
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/bindings/
  )

if (${ENABLE_DOC})
  add_subdirectory(
    ${CMAKE_SOURCE_DIR}/docs/
    )
endif (${ENABLE_DOC})

# compile executables
file( GLOB binaries "${CMAKE_SOURCE_DIR}/examples/*.cc")
foreach(binaryfile ${binaries})
  get_filename_component(binaryname  ${binaryfile} NAME_WE)
  add_executable(${binaryname} ${binaryfile})
  target_link_libraries(${binaryname} ${Boost_LIBRARIES} rascal)
endforeach(binaryfile ${binaries})

# or copy them
file (GLOB pybins "${CMAKE_SOURCE_DIR}/bin/*.py")
foreach(pybin ${pybins})
  get_filename_component(binaryname  ${pybin} NAME_WE)
  configure_file(
    ${pybin}
    "${CMAKE_BINARY_DIR}/${binaryname}.py"
    COPYONLY)
endforeach(pybin ${pybins})

# Checking hh/cc files for compliance with coding conventions
if(CPPLINT_FOUND)
  # register the folders to check for cpp code conventions compliance
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/src" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/tests" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/bindings" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/examples"
  "--filter=-build/namespaces")
endif()

# Add all hh/cc files for autoformatting with clang-format
if(CLANG_FORMAT_FOUND)
  # register the folders to apply clang-format with given configuration
  # .clang-format in project root folder
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/src" "")
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/tests" "")
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/bindings" "")
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/examples")
endif()

# Add python files for autoformatting with autopep8
if(AUTOPEP8_FOUND)
  # register the folders to apply clang-format with given configuration
  # .clang-format in project root folder
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/src" "")
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/tests" "")
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/bindings" "")
  clang_format_add_subdirectory("${CMAKE_SOURCE_DIR}/examples")
endif()

################################################################################
# make python library of rascal usable in build folder
add_custom_target(pyrascalcopy ALL)
add_custom_command(TARGET pyrascalcopy POST_BUILD
  COMMAND cp -R
  ${CMAKE_SOURCE_DIR}/rascal ${CMAKE_BINARY_DIR}
  )
