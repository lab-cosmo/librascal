cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

#set(LIBTORCH_PATH "${LIBTORCH_PATH}")
#message(FATAL_ERROR "LIBTORCH_PATH  ${LIBTORCH_PATH}")
set(Torch_DIR "${LIBTORCH_PATH}/share/cmake/Torch/")
# not sure if this line works properly
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH} ${LIBTORCH_PATH}")
find_package(Torch REQUIRED)
# not sure if this line is needed
include_directories(${TORCH_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}" CACHE INTERNAL "")

# ENV is used because standard variables are not cached outsite of this cmake file
set(ENV{TORCH_LIBRARIES} "${TORCH_LIBRARIES}")
set(ENV{LIBTORCH_PATH} "${LIBTORCH_PATH}")
#message(FATAL_ERROR "TORCH_LIBRARIES $ENV{TORCH_LIBRARIES}")
#message(FATAL_ERROR "LIBTORCH_PATH ${LIBTORCH_PATH}")
