#
# CMake module to C++ static analysis against
# Google C++ Style Guide (https://google.github.io/styleguide/cppguide.html)
#
# For more detials please follow links:
#
# - https://github.com/google/styleguide
# - https://pypi.python.org/pypi/cpplint
# - https://github.com/theandrewdavis/cpplint
#
# Copyright (c) 2016 Piotr L. Figlarek
#
# Usage ----- Include this module via CMake include(...) command and then add
# each source directory via introduced by this module
# cpplint_add_subdirectory(...) function. Added directory will be recursivelly
# scanned and all available files will be checked.
#
# Example
# -------
# # include CMake module
# include(cmake/cpplint.cmake)
#
# # add all source code directories
# cpplint_add_subdirectory(core)
# cpplint_add_subdirectory(modules/c-bind)
#
# License (MIT)
# -------------
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Custom target to invoke clang-format on all c++ files in the project
# The style is customized in the /.clang-format file in project root

# target to run cpplint.py for all configured sources
set(CLANG_FORMAT_TARGET pretty-cpp)

# project root directory
set(CLANG_FORMAT_PROJECT_ROOT
  ${PROJECT_SOURCE_DIR} CACHE STRING "Project ROOT directory"
  )

mark_as_advanced(CLANG_FORMAT_PROJECT_ROOT)
mark_as_advanced(CLANG_FORMAT)

# find clang-format executable
find_program(CLANG_FORMAT name clang-format HINTS $ENV/usr/bin)
if(CLANG_FORMAT)
    message(STATUS "clang-format executable: ${CLANG_FORMAT}")
    # common target to concatenate all cpplint.py targets
    add_custom_target(${CLANG_FORMAT_TARGET})
    set(CLANG_FORMAT_FOUND TRUE)
else()
    message(STATUS "The optional clang-format auto formatter has not been"
                   " found. For more information see"
                   " https://clang.llvm.org/docs/ClangFormat.html")
    set(CLANG_FORMAT_FOUND FALSE)
    return()
endif()

execute_process(
    COMMAND ${CLANG_FORMAT} --version
    OUTPUT_VARIABLE CLANG_FORMAT_VERSION
)
string(REGEX REPLACE ".*clang-format version ([0-9.]+).*" "\\1" CLANG_FORMAT_VERSION "${CLANG_FORMAT_VERSION}")
if("${CLANG_FORMAT_VERSION}" VERSION_LESS "8.0")
    message(STATUS "Unsupported clang-format (version ${CLANG_FORMAT_VERSION}): "
                   "a more recent version (at least 8.0) is required")
    set(CLANG_FORMAT_FOUND FALSE)
    return()
endif()

# use clang-format to autoformat source code files DIR
function(clang_format_add_subdirectory DIR)
    # create relative path to the directory
    set(ABSOLUTE_DIR ${DIR})

    set(EXTENSIONS cc,hh)
    set(FILES_TO_CHECK ${FILES_TO_CHECK}
      ${ABSOLUTE_DIR}/*.cc ${ABSOLUTE_DIR}/*.hh
      )

    # find all source files inside project
    file(GLOB_RECURSE LIST_OF_FILES ${FILES_TO_CHECK})

    # create valid target name for this check
    string(REGEX REPLACE "/" "." TEST_NAME ${DIR})
    set(TARGET_NAME ${CLANG_FORMAT_TARGET}.${TEST_NAME})

    # autoformat all hh/cc files in project
    add_custom_target(${TARGET_NAME}
        COMMAND ${CLANG_FORMAT} -i ${LIST_OF_FILES}
        DEPENDS ${LIST_OF_FILES}
        COMMENT "clang-format: autoformatting all hh/cc files"
    )

    # run this target when `make pretty` is invoked
    add_dependencies(${CLANG_FORMAT_TARGET} ${TARGET_NAME})

endfunction()
