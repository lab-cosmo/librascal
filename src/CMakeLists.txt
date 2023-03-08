set(RASCAL_SOURCES
    rascal/utils/json_io.cc
    rascal/utils/units.cc
    rascal/utils/utils.cc

    rascal/math/bessel.cc
    rascal/math/hyp1f1.cc
    rascal/math/interpolator.cc
    rascal/math/gauss_legendre.cc
    rascal/math/spherical_harmonics.cc
    rascal/math/kvec_generator.cc
    rascal/structure_managers/structure_manager_lammps.cc
    rascal/structure_managers/structure_manager_centers.cc
    rascal/representations/calculator_base.cc
)

add_library(${LIBRASCAL_NAME} ${RASCAL_SOURCES})

target_include_directories(${LIBRASCAL_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${LIBRASCAL_NAME} PUBLIC Eigen3::Eigen)
target_link_libraries(${LIBRASCAL_NAME} PUBLIC ${WIGXJPF_NAME})

install(TARGETS ${LIBRASCAL_NAME} DESTINATION lib)
