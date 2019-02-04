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

# target to run autopep8.py for all configured sources
set(AUTOPEP8_TARGET prettypython CACHE STRING "Name of python autoformatter")

# project root directory
set(AUTOPEP8_PROJECT_ROOT ${PROJECT_SOURCE_DIR}
  CACHE STRING "Project ROOT directory"
  )

# find autopep8.py script
find_file(AUTOPEP8 name autopep8 HINTS $ENV{HOME}.local/bin)
if(AUTOPEP8)
    message(STATUS "autopep8 auto formatter: ${AUTOPEP8}")
    # common target to concatenate all autopep8.py targets
    add_custom_target(${AUTOPEP8_TARGET})
    set(AUTOPEP8_FOUND TRUE)
else()
    message(STATUS "The optional autopep8 parser has not been found. "
      "For more information see"
      " https://github.com/hhatto/autopep8")
    set(AUTOPEP8_FOUND FALSE)
endif()




# use autopep8.py to check source code files inside DIR directory
function(autopep8_add_subdirectory DIR)
    # create relative path to the directory
    set(ABSOLUTE_DIR ${DIR})

    set(EXTENSIONS py)
    set(FILES_TO_CHECK ${FILES_TO_CHECK}
      ${ABSOLUTE_DIR}/*.py
      )

    # find all source files inside project
    file(GLOB_RECURSE LIST_OF_FILES ${FILES_TO_CHECK})

    # create valid target name for this check
    string(REGEX REPLACE "/" "." TEST_NAME ${DIR})
    set(TARGET_NAME ${AUTOPEP8_TARGET}.${TEST_NAME})

    # perform autopep8 check
    add_custom_target(${TARGET_NAME}
        COMMAND ${AUTOPEP8} "--in-place"
                           ${LIST_OF_FILES}
        DEPENDS ${LIST_OF_FILES}
        COMMENT "autopep8: Checking source code style"
    )

    # run this target, when `prettypython` is
    add_dependencies(${AUTOPEP8_TARGET} ${TARGET_NAME})

endfunction()
