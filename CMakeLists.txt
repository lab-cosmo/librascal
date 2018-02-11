cmake_minimum_required(VERSION 2.8.12)

project(proteus)

option(NDARRAY_EIGEN "Enable Eigen" ON)

# put our local cmake find scripts at the beginning of the cmake
# module search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# enable C++14 support
add_definitions(-std=c++14)

### Eigen dependency tests
if(NDARRAY_EIGEN)
    find_package(Eigen3 REQUIRED)
    include_directories( ${EIGEN3_INCLUDE_DIR} )
endif(NDARRAY_EIGEN)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_compile_options(-msse2 )

add_subdirectory(pybind11)


set(PROTEUS_SRC_DIR ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(${PROTEUS_SRC_DIR})



#pybind11_add_module(proteus src/module.cc)
#add_subdirectory(proteus)