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


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
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

# Try to find optional cpp linter
include(cpplint)

include_directories(
  ${CMAKE_SOURCE_DIR}/tests
  ${CMAKE_SOURCE_DIR}/src
 )

include_directories(SYSTEM
  ${CMAKE_SOURCE_DIR}/external
  )

set(ENABLE_DOC "OFF" CACHE BOOL
  "If on, a target 'dev_doc' (which builds the documentation) is added")
set(BUILD_TESTS "OFF" CACHE BOOL
  "If on, builds the unit tests")
set(SUPPRESS_LIBRARY_WARNINGS "ON" CACHE BOOL
  "Suppress warnings from  Boost library")

set(INSTALL_PATH "" CACHE STRING
  "Path to install the python library. By default installs in the detected python lib directory.")

set(USER_INSTALL "OFF" CACHE BOOL
  "Switch installation from the main python site-packages to the user site-packages")


  find_package(PythonLibsNew 3 MODULE REQUIRED)

  # Get where python site-packages is
  execute_process ( COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import \
                    get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # Get where python local site-packages is
  execute_process ( COMMAND "${PYTHON_EXECUTABLE}" -c "import site; \
                     print(site.USER_SITE)"
    OUTPUT_VARIABLE PYTHON_USER_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # Default install in site package of the detected python library
  if("${INSTALL_PATH}" STREQUAL "" AND NOT ${USER_INSTALL})
    set(CMAKE_INSTALL_PREFIX "${PYTHON_SITE_PACKAGES}/rascal/"
      CACHE PATH "..." FORCE)

  # install in the user's site package of the detected python library
  elseif("${INSTALL_PATH}" STREQUAL "" AND ${USER_INSTALL})
    set(CMAKE_INSTALL_PREFIX "${PYTHON_USER_SITE_PACKAGES}/rascal/"
    CACHE PATH "..." FORCE)

  # install in the specified path
  elseif(NOT "${INSTALL_PATH}" STREQUAL "")
    set(CMAKE_INSTALL_PREFIX "${INSTALL_PATH}/rascal/" CACHE PATH "..." FORCE)

  endif()

  set(CMAKE_SKIP_BUILD_RPATH  FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif()

  message(STATUS "######### ${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")

add_subdirectory(
    ${CMAKE_SOURCE_DIR}/bindings/
    )

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
  add_test(python_binding_test python_binding_tests.py)

  #build tests
  file( GLOB TEST_SRCS "${CMAKE_SOURCE_DIR}/tests/test_*.cc")

  add_executable(main_test_suite tests/main_test_suite.cc ${TEST_SRCS})

  # https://stackoverflow.com/questions/7781374/boost-test-linking
  ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

  target_include_directories(main_test_suite SYSTEM PRIVATE ${Boost_INCLUDE_DIR})
  target_link_libraries(main_test_suite ${Boost_LIBRARIES} rascal)
  add_test(main_test_suite main_test_suite --report_level=detailed --build_info=TRUE)
  configure_file("tests/alanine-X.json" ${CMAKE_BINARY_DIR} COPYONLY)
  configure_file("tests/simple_cubic_9.json" ${CMAKE_BINARY_DIR} COPYONLY)

  file(COPY "tests/reference_data" DESTINATION ${CMAKE_BINARY_DIR})

endif(${BUILD_TESTS})

# compile the library
# add_compile_options(-Werror) # warnings are errors!
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/src/
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


if(CPPLINT_FOUND)
  # register the folders to check for cpp code conventions compliance
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/src" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/tests" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/bindings" "")
  cpplint_add_subdirectory("${CMAKE_SOURCE_DIR}/examples"
  "--filter=-build/namespaces")
endif()
