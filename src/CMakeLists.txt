


set(proteus_SRC
    module.cc)


add_library(proteus ${proteus_SRC})
target_link_libraries(proteus Eigen3::Eigen)
