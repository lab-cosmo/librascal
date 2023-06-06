# Downloads and compiles the wigxjpf library staticaly and link it into rascal

set(SOURCE_DIR "${_wigxjpf_external_dir}/${package}")

download_external_project(wigxjpf
    URL "https://github.com/lab-cosmo/wigxjpf/archive/refs/tags/${_wigxjpf_version}.tar.gz"
    BACKEND NONE
    THIRD_PARTY_SRC_DIR ${_wigxjpf_external_dir}
    BUILD_COMMAND "true"
    ${_wigxjpf_update}
)

file(WRITE ${SOURCE_DIR}/cfg/wigxjpf_auto_config.h "")

set(WIGXJPF_SOURCES
    ${SOURCE_DIR}/src/c_wrap.c
    ${SOURCE_DIR}/src/calc.c
    ${SOURCE_DIR}/src/prime_factor.c
    ${SOURCE_DIR}/src/trivial_zero.c
)
add_library("${WIGXJPF_NAME}" OBJECT ${WIGXJPF_SOURCES})

target_include_directories("${WIGXJPF_NAME}" PUBLIC
    ${SOURCE_DIR}/inc
    ${SOURCE_DIR}/cfg
)

set_target_properties("${WIGXJPF_NAME}" PROPERTIES POSITION_INDEPENDENT_CODE ON)
