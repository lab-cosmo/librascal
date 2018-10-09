add_library(rascal-cephes "")

add_subdirectory(cephes)

install(TARGETS rascal-cephes DESTINATION lib)