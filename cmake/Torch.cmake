cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(LIBTORCH_PATH "/ssd/local/code/libtorch/" )
#set(TORCH_LIBRARIES "/ssd/local/code/libtorch/lib")
set(Torch_DIR "${LIBTORCH_PATH}/share/cmake/Torch/")
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH} ${LIBTORCH_PATH}")
set(TORCH_INCLUDE_DIRS "${LIBTORCH_PATH}/include")
#set(Torch_DIR ${LIBTORCH_PATH})
find_package(Torch REQUIRED)
include_directories(${TORCH_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
#add_executable(torch_example "${PROJECT_SOURCE_DIR}/examples/cpp/torch_example.cc")
#target_link_libraries(torch_example "${TORCH_LIBRARIES}")
#set_property(TARGET torch_example PROPERTY CXX_STANDARD 14)

#find_package(Torch REQUIRED)
#add_subdirectory("${PROJECT_SOURCE_DIR}/pytorch")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

#add_subdirectory(${_Torch_external_dir}/build)


#add_executable(example-app example-app.cpp)
#target_link_libraries(example-app "${TORCH_LIBRARIES}")
#set_property(TARGET example-app PROPERTY CXX_STANDARD 14)
