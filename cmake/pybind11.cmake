download_external_project(pybind11
  URL "https://github.com/pybind/pybind11.git"
  TAG "v${_pybind11_version}"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${_pybind11_external_dir}
  NO_UPDATE
  )

set(PYBIND11_LTO_CXX_FLAGS "" CACHE INTERNAL "")
set(PYBIND11_LTO_LINKER_FLAGS "" CACHE INTERNAL "")

add_subdirectory(${_pybind11_external_dir}/pybind11)

include_directories(SYSTEM ${PYBIND11_INCLUDE_DIR})

mark_as_advanced_prefix(PYBIND11)
mark_as_advanced(USE_PYTHON_INCLUDE_DIR)

set(${package}_FOUND TRUE CACHE INTERNAL "To avoid cyclic search" FORCE)
set(${package}_FOUND_EXTERNAL TRUE CACHE INTERNAL "" FORCE)
