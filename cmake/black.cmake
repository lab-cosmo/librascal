
set(BLACK_TARGET pretty-python)
set(BLACK_LINT_TARGET pylint)

# project root directory
set(BLACK_PROJECT_ROOT ${PROJECT_SOURCE_DIR}
  CACHE STRING "Project ROOT directory"
  )

mark_as_advanced(BLACK_PROJECT_ROOT)
mark_as_advanced(BLACK)

# find black.py script
find_file(BLACK name black HINTS $ENV{HOME}.local/bin)
if(BLACK)
    message(STATUS "BLACK auto formatter: ${BLACK}")
    # common target to concatenate all BLACK.py targets
    add_custom_target(${BLACK_TARGET})
    add_custom_target(${BLACK_LINT_TARGET})
    set(BLACK_FOUND TRUE)

    # run black on python files that are at the root directory
    set(ROOT_DIR_FILES "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
    add_custom_target(black.setup
        COMMAND ${BLACK} "-l 88" "-v" ${ROOT_DIR_FILES}
        COMMENT "black: reformatting python source code"
    )
    add_custom_target(black.lint.setup
        COMMAND ${BLACK} "--check" "-l 88" "-v" ${ROOT_DIR_FILES}
        COMMENT "black: Checking python source code style"
    )
    add_dependencies(${BLACK_TARGET} black.setup)
    add_dependencies(${BLACK_LINT_TARGET} black.lint.setup)


else()
    message(STATUS "The optional black formater has not been found. "
      "For more information see"
      " https://github.com/psf/black")
    set(BLACK_FOUND FALSE)
endif()

# use black.py to check source code files inside DIR directory
function(black_add_subdirectory DIR)
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
    set(TARGET_NAME ${BLACK_TARGET}.${TEST_NAME})
    set(LINT_TARGET_NAME lint.${BLACK_TARGET}.${TEST_NAME})

    # perform black reformatting
    add_custom_target(${TARGET_NAME}
        COMMAND ${BLACK} "-l 88" "-v" ${LIST_OF_FILES}
        DEPENDS ${LIST_OF_FILES}
        COMMENT "black: reformatting python source code"
    )

    # perform black reformatting
    add_custom_target(${LINT_TARGET_NAME}
        COMMAND ${BLACK} "--check" "-l 88" "-v" ${LIST_OF_FILES}
        DEPENDS ${LIST_OF_FILES}
        COMMENT "black: Checking python source code style"
    )

    # run this target, when `pretty-python` is
    add_dependencies(${BLACK_TARGET} ${TARGET_NAME})
    # run this target, when `pylint` is
    add_dependencies(${BLACK_LINT_TARGET} ${LINT_TARGET_NAME})

endfunction()