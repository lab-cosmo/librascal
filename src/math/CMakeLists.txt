add_library(cephes "")

add_subdirectory(cephes)

install(TARGETS cephes DESTINATION lib)



add_library(rascal-math "")

install(TARGETS rascal-math DESTINATION lib)

set (math_SRC
  math_utils.cc
  math_interface.cc
  )

target_sources(rascal-math PRIVATE ${math_SRC})

target_link_libraries(rascal-math PRIVATE cephes)