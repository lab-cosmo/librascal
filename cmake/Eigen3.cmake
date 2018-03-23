download_external_project(Eigen3
  URL "https://github.com/eigenteam/eigen-git-mirror.git"
  TAG "${_Eigen3_version}"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${_Eigen3_external_dir}
  ${_Eigen3_update}
  )

file(READ "${_Eigen3_external_dir}/Eigen3/Eigen/src/Core/util/Macros.h" _eigen_version_header)
string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen_world_version_match "${_eigen_version_header}")
set(_eigen_world_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen_major_version_match "${_eigen_version_header}")
set(_eigen_major_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen_minor_version_match "${_eigen_version_header}")
set(_eigen_minor_version "${CMAKE_MATCH_1}")
set(Eigen3_VERSION ${_eigen_world_version}.${_eigen_major_version}.${_eigen_minor_version})

set(EIGEN3_INCLUDE_DIR ${_Eigen3_external_dir}/Eigen3 CACHE PATH "Eigen include directory")
add_library(eigen3 INTERFACE)
add_library(Eigen3::Eigen ALIAS eigen3)
target_include_directories(eigen3 SYSTEM INTERFACE
  $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>)

set(Eigen3_FOUND TRUE CACHE BOOL INTERNAL "")
set(Eigen3_FOUND_EXTERNAL TRUE CACHE BOOL INTERNAL "")
message(STATUS "Eigen3 ${Eigen3_VERSION}")
