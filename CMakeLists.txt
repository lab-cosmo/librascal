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
# Copyright © 2018 Till Junge, Felix Musil, Markus Stricker, COSMO (EPFL), LAMMM
# (EPFL)
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

# https://cmake.org/pipermail/cmake/2008-September/023808.html
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their
# definition
# and dump it in the cache along with proper documentation, otherwise set
# CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()
# Set the possible values of build type for cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
"Release" "Debug" "RelWithDebInfo" "MinSizeRel")

project(rascal)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS ON)
set(RASCAL_PYTHON_MAJOR_VERSION 3)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(rascalTools)

########## MAIN VARIABLE TO CHANGE THE BUILD AND INSTALL OF THE LIBRARY ########
set(BUILD_BINDINGS "ON" CACHE BOOL
  "If on, the python bindings will be built.")

set(USER "OFF" CACHE BOOL
  "Switch installation folder")

set(INSTALL_PATH "" CACHE STRING
  "Path to install the libraries")

set(ENABLE_DOC "OFF" CACHE BOOL
  "If on, a target 'dev_doc' (which builds the documentation) is added.")

set(BUILD_TESTS "OFF" CACHE BOOL
  "If on, builds the unit tests")

set(ENABLE_PROFILING "OFF" CACHE BOOL
  "If on, adds gprof profiling information to some of the examples. \
  (GNU compilers only)")

message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")

########## COMPILATION FLAGS ##########
add_compile_options(-Wall -Wextra)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
# to avoid conflicts with built-in functions
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-builtin")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
    OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  add_compile_options(-Wno-missing-braces)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
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

########## START CONFIGURATION ##########

# set the name of the c++ libraries
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR
    CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(LIBRASCAL_NAME "rascal-dbg")
  set(LIBCEPHES_NAME "cephes-dbg")
else()
  set(LIBRASCAL_NAME "rascal")
  set(LIBCEPHES_NAME "cephes")
endif()

add_external_package(Eigen3 VERSION 3.3.4 CONFIG)

# Try to find optional cpp linter
include(cpplint)

# Check for clang-format
include(clangformat)

# CHeck for autopep8
include(autopep8)

include_directories(
  ${CMAKE_SOURCE_DIR}/tests
  ${CMAKE_SOURCE_DIR}/src
 )

include_directories(SYSTEM
  ${CMAKE_SOURCE_DIR}/external
  )


# set up bindings
if(BUILD_BINDINGS)
  find_package(PythonLibsNew 3 MODULE REQUIRED)

  # warning about potential side effects of the installation of the python library
  if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR
      CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    if(NOT PYTHON_IS_DEBUG)
      message(WARNING "The current configuration will build debug libraries "
              "not linked to a debug python interpreter. "
              "Installation target might overwrite the python bindings "
              "library with a debug build !"
              "see https://wiki.ubuntu.com/PyDbgBuilds or "
              "https://gist.github.com/bombs-kim/c9848c2b09962f2fd753b48b6d2cd87ffor sources and "
              "https://jml.io/2015/08/debugging-python-with-gdb.html for more informations "
              " on how to debug c++ code from python.")
    endif()
  endif()

  add_external_package(pybind11 VERSION 2.2.2 CONFIG)

    # Get where python site-packages is
    if(NOT USER)
    execute_process ( COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import \
                      get_python_lib; print(get_python_lib())"
      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  else()
    # Get where python local site-packages is
    execute_process ( COMMAND "${PYTHON_EXECUTABLE}" -c "import site; \
                      print(site.USER_SITE)"
      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  endif()

  execute_process ( COMMAND "${PYTHON_EXECUTABLE}" -c "import sys; \
                    print(sys.prefix)"
    OUTPUT_VARIABLE PACKAGE_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

# figure out where to install the library
if(NOT "${INSTALL_PATH}" STREQUAL "")
  set(CMAKE_INSTALL_PREFIX "${INSTALL_PATH}")
elseif(NOT BUILD_BINDINGS AND NOT USER)
  set(CMAKE_INSTALL_PREFIX "/usr")
  message(STATUS "INSTALLATION REQUIRES SUDO RIGHTS")
elseif(NOT BUILD_BINDINGS AND USER)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local")
elseif(BUILD_BINDINGS AND NOT USER)
  set(CMAKE_INSTALL_PREFIX "${PACKAGE_PREFIX}")
elseif(BUILD_BINDINGS AND USER)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local")
endif()

message(STATUS "Installation ROOT: ${CMAKE_INSTALL_PREFIX}")


# Always use full RPATH
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# utility to force copy of files...
add_custom_command(
    OUTPUT always_rebuild
    COMMAND ${CMAKE_COMMAND} -E echo
  )

# set up the tests
if (${BUILD_TESTS})
  find_package(Boost REQUIRED COMPONENTS unit_test_framework )
  enable_testing()
  # build c++ tests and add to test suite
  file( GLOB TEST_SRCS "${CMAKE_SOURCE_DIR}/tests/test_*.cc")

  add_executable(main_test_suite tests/main_test_suite.cc ${TEST_SRCS})

  # https://stackoverflow.com/questions/7781374/boost-test-linking
  ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

  target_include_directories(main_test_suite SYSTEM PRIVATE ${Boost_INCLUDE_DIR})
  target_link_libraries(main_test_suite ${Boost_LIBRARIES} "${LIBRASCAL_NAME}")
  add_test(main_test_suite main_test_suite --report_level=detailed --build_info=TRUE)
  configure_file("tests/alanine-X.json" ${CMAKE_BINARY_DIR} COPYONLY)
  configure_file("tests/simple_cubic_9.json" ${CMAKE_BINARY_DIR} COPYONLY)


  if(BUILD_BINDINGS)

    add_custom_target(py_test ALL
      DEPENDS
        copy_test_py
      )

    add_dependencies(py_test main_test_suite)

    add_custom_command(OUTPUT copy_test_py
      COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_SOURCE_DIR}/tests/python ${CMAKE_BINARY_DIR}/
      DEPENDS
          always_rebuild
    )

    add_test(python_binding_test python_binding_tests.py)
  endif()


endif(${BUILD_TESTS})

# copy reference data that might be used in the example
add_custom_target(copy_reference_data ALL
    DEPENDS
    copy_test_reference_data
  )

add_custom_command(OUTPUT copy_test_reference_data
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/tests/reference_data ${CMAKE_BINARY_DIR}/reference_data
    DEPENDS
        always_rebuild
  )

# put it there to avoid warnings from boost tests
# add_compile_options(-Werror) # warnings are errors!

# compile the library
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/src/
  )

# set up bindings
if(BUILD_BINDINGS)
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/bindings/
)
endif()

# set up documentation
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
  target_link_libraries(${binaryname} "${LIBRASCAL_NAME}")
endforeach(binaryfile ${binaries})
configure_file("${CMAKE_SOURCE_DIR}/examples/crystal_structure.json" ${CMAKE_BINARY_DIR} COPYONLY)
configure_file("${CMAKE_SOURCE_DIR}/examples/alanine-X.json" ${CMAKE_BINARY_DIR} COPYONLY)

# Add profiling information for the appropriate examples
if (ENABLE_PROFILING)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options("soap_profile" PRIVATE -ggdb)
    # target_compile_options("soap_profile" PRIVATE -pg)
    # set_target_properties("soap_profile" PROPERTIES LINK_FLAGS -pg)
  endif()
endif()



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
  # register the folders to apply autopep8 with given configuration
  # .pycodestyle in project root folder
  autopep8_add_subdirectory("${CMAKE_SOURCE_DIR}/tests" "")
  autopep8_add_subdirectory("${CMAKE_SOURCE_DIR}/bindings" "")
endif()
