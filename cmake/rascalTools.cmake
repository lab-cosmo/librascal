#==============================================================================
# file   RascalTools
#
# @author Nicolas Richart <nicolas.,richart@epfl.ch>
# @author Till Junge <till.junge@epfl.ch>
#
# @date   11 Jan 2018
#
# @brief  some tool to help to do stuff with cmake in rascal
#
# @section LICENSE
#
# Copyright  2018 Till Junge
#
# rascal is free software; HANDLE LICENSES HERE
#

#[=[.rst:
rascalTools
-------------

This module provide some helper functions for rascal

::

download_external_project(project_name
  URL <url>
  BACKEND <GIT|HG|NONE>
  [TAG <tag>]
  [THIRD_PARTY_SRC_DIR <dir>]
  [NO_UPDATE]
  )

Downloads the external project based on the uri

``THIRD_PARTY_SRC_DIR <dir>``
  Specifies where to download the source

``NO_UPDATE``
  Does not try to update existing download

::

mark_as_advanced_prefix(prefix)

Marks as advanced all variables whoes names start with prefix

::

add_external_dep(package
  [IGNORE_SYSTEM]
  [VERSION <version>]
  [...])

Tries to find the external package and if not found but a local
${package}.cmake files exists, it includes it. The extra arguments are passed
to find_package()

``IGNORE_SYSTEM`` does not do the find_package

``VERSION`` Specifies the required minimum version

]=]
#

function(download_external_project project_name)
  include(CMakeParseArguments)

  set(_dep_flags
    NO_UPDATE)
  set(_dep_one_variables
    URL
    TAG
    BACKEND
    THIRD_PARTY_SRC_DIR
    CONFIGURE_COMMAND
    BUILD_COMMAND
    INSTALL_COMMAND
    )
  set(_dep_multi_variables)
  cmake_parse_arguments(_dep_args
    "${_dep_flags}"
    "${_dep_one_variables}"
    "${_dep_multi_variables}"
    ${ARGN}
    )

  if(NOT _dep_args_URL)
    message(FATAL_ERROR "You have to provide a URL for the project ${project_name}")
  endif()

  if(NOT _dep_args_BACKEND)
    message(FATAL_ERROR "You have to provide a backend to download ${project_name}")
  endif()


  set(_src_dir ${PROJECT_SOURCE_DIR}/third-party/${project_name})
  if (_dep_args_THIRD_PARTY_SRC_DIR)
    set(_src_dir ${_dep_args_THIRD_PARTY_SRC_DIR}/${project_name})
  endif()

  if(EXISTS ${_src_dir} AND _dep_args_NO_UPDATE)
    return()
  endif()

  set(_working_dir ${PROJECT_BINARY_DIR}/third-party/${project_name}-download)

  if(_dep_args_BACKEND STREQUAL "NONE")
    file(WRITE ${_working_dir}/CMakeLists.txt
  "
cmake_minimum_required(VERSION 3.1)
project(${project_name}-download NONE)
include(ExternalProject)
ExternalProject_Add(${project_name}
SOURCE_DIR ${_src_dir}
BINARY_DIR ${_working_dir}
URL ${_dep_args_URL}
CONFIGURE_COMMAND \"\"
BUILD_COMMAND ${_dep_args_BUILD_COMMAND}
INSTALL_COMMAND   \"\"
TEST_COMMAND      \"\"
)
")
#
  else()
    if(_dep_args_TAG)
      set(_ep_tag "${_dep_args_BACKEND}_TAG ${_dep_args_TAG}")
    endif()

    file(WRITE ${_working_dir}/CMakeLists.txt
        "
cmake_minimum_required(VERSION 3.1)
project(${project_name}-download NONE)
include(ExternalProject)
ExternalProject_Add(${project_name}
  SOURCE_DIR ${_src_dir}
  BINARY_DIR ${_working_dir}
  ${_dep_args_BACKEND}_REPOSITORY ${_dep_args_URL}
  ${_ep_tag}
  CONFIGURE_COMMAND \"\"
  BUILD_COMMAND     \"\"
  INSTALL_COMMAND   \"\"
  TEST_COMMAND      \"\"
  )
")

  endif(_dep_args_BACKEND STREQUAL "NONE")


  message(STATUS "Downloading ${project_name} ${_dep_args_GIT_TAG}")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE _result
    WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/configure-out.log
    ERROR_FILE ${_working_dir}/configure-error.log)

  if(_result)
    message(FATAL_ERROR "Something went wrong (${_result}) during the download"
      " process of ${project_name} check the file"
      " ${_working_dir}/configure-error.log for more details")
  endif()

  execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    RESULT_VARIABLE _result
    WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/build-out.log
    ERROR_FILE ${_working_dir}/build-error.log)

  if(_result)
    message(FATAL_ERROR "Something went wrong (${_result}) during the download"
      " process of ${project_name} check the file"
      " ${_working_dir}/build-error.log for more details")
  endif()

endfunction()

# ------------------------------------------------------------------------------
function(mark_as_advanced_prefix prefix)
  get_property(_list DIRECTORY PROPERTY VARIABLES)
  foreach(_var ${_list})
    if(${_var} MATCHES "^${prefix}.*")
      mark_as_advanced(${_var})
    endif()
  endforeach()
endfunction()

# ------------------------------------------------------------------------------
function(add_external_package package)
  include(CMakeParseArguments)

  set(_cmake_includes ${PROJECT_SOURCE_DIR}/cmake)
  set(_${package}_external_dir ${PROJECT_BINARY_DIR}/external)

  set(_aep_flags
    IGNORE_SYSTEM
    )
  set(_aep_one_variables
    VERSION
    )
  set(_aep_multi_variables)

  cmake_parse_arguments(_aep_args
    "${_aep_flags}"
    "${_aep_one_variables}"
    "${_aep_multi_variables}"
    ${ARGN}
    )

  if(_aep_args_VERSION)
    set(_${package}_version ${_aep_args_VERSION})
  endif()

  if(NOT EXISTS ${_cmake_includes}/${package}.cmake)
    set(_required REQUIRED)
  endif()

  if(NOT _aep_args_IGNORE_SYSTEM)
    find_package(${package} ${_${package}_version} ${_required} ${_aep_UNPARSED_ARGUMENTS} QUIET)
    if(${package}_FOUND AND NOT ${package}_FOUND_EXTERNAL)
      return()
    endif()
  endif()



  if(EXISTS ${_cmake_includes}/${package}.cmake)
    include(${_cmake_includes}/${package}.cmake)
  endif()
endfunction()
