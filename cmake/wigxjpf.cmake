# Downloads and compiles the wigxjpf library staticaly and link it into rascal

set(SOURCE_DIR "${_wigxjpf_external_dir}/${package}")

download_external_project(wigxjpf
    URL "http://fy.chalmers.se/subatom/wigxjpf/wigxjpf-${_wigxjpf_version}.tar.gz"
    BACKEND NONE
    THIRD_PARTY_SRC_DIR ${_wigxjpf_external_dir}
    BUILD_COMMAND "make -C ${SOURCE_DIR}"
    ${_wigxjpf_update}
)

add_library("${WIGXJPF_NAME}" STATIC IMPORTED)

set_target_properties("${WIGXJPF_NAME}" PROPERTIES
    IMPORTED_LOCATION ${SOURCE_DIR}/lib/libwigxjpf.a
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/inc;${SOURCE_DIR}/cfg"
)

add_dependencies("${WIGXJPF_NAME}" wigxjpf)
