download_external_project(boost_test
  URL "https://github.com/boostorg/test.git"
  TAG "boost-${_boost_test_version}"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${_boost_test_external_dir}
  NO_UPDATE
  )


set(Boost_NO_SYSTEM_PATHS ON)
set(BOOST_ROOT ${_boost_test_external_dir})
set(boost_test_INCLUDE_DIR ${_boost_test_external_dir}/include CACHE PATH "boost test include directory")

include_directories(${boost_test_INCLUDE_DIR})

find_package(Boost ${_boost_test_version} COMPONENTS unit_test_framework)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  #add_executable(foo foo.cc)
endif()

add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

message(STATUS "Boost Test" ${Boost_FOUND})

message(STATUS "Boost Test" ${BOOST_SOURCE})
message(STATUS "Boost Test" ${Boost_LIBRARY_DIRS})

