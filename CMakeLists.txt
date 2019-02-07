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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS ON)
set(RASCAL_PYTHON_MAJOR_VERSION 3)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(rascalTools)

########## MAIN VARIABLE TO CHANGE THE BUILD AND INSTALL OF THE LIBRARY ############
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

message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")

########## COMPILATION FLAGS ##########
add_compile_options(-Wall -Wextra)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
# to avoid conflicts with built-in functions
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-builtin")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
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
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(LIBRASCAL_NAME "rascal-dbg")
  set(LIBCEPHES_NAME "cephes-dbg")
else()
  set(LIBRASCAL_NAME "rascal")
  set(LIBCEPHES_NAME "cephes")
endif()

add_external_package(Eigen3 VERSION 3.3.4 CONFIG)

# Try to find optional cpp linter
include(cpplint)

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
  if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
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


# set up the tests
if (${BUILD_TESTS})
  find_package(Boost REQUIRED COMPONENTS unit_test_framework )
  enable_testing()
  # copy python test

  file( GLOB PY_TEST_SRCS "${CMAKE_SOURCE_DIR}/tests/python_*.py")

  foreach(pytest ${PY_TEST_SRCS})
    get_filename_component(pytest_name ${pytest} NAME)
    configure_file(
      ${pytest}
      "${CMAKE_BINARY_DIR}/${pytest_name}"
      COPYONLY)
  endforeach(pytest ${PY_TEST_SRCS})

  if(BUILD_BINDINGS)
    add_test(python_binding_test python_binding_tests.py)
  endif()

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
  file(COPY "tests/reference_data" DESTINATION ${CMAKE_BINARY_DIR})

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



# or copy them
file (GLOB pybins "${CMAKE_SOURCE_DIR}/bin/*.py")
foreach(pybin ${pybins})
  get_filename_component(binaryname  ${pybin} NAME_WE)
  configure_file(
    ${pybin}
    "${CMAKE_BINARY_DIR}/${binaryname}.py"
    COPYONLY)
endforeach(pybin ${pybins})


if(CPPLINT_FOUND)
  # register the folders to check for cpp code conventions compliance
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/src" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/tests" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/bindings" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/examples"
  "--filter=-build/namespaces")
endif()
