add_executable(kappa_a kappa_a.cpp)
target_link_libraries(kappa_a mudiraclib)

add_executable(calc_state calc_state.cpp)
target_link_libraries(calc_state mudiraclib)

add_executable(size_effect size_effect.cpp)
target_link_libraries(size_effect mudiraclib)

add_executable(eigenenergies eigenenergies.cpp)
target_link_libraries(eigenenergies mudiraclib)

add_executable(hydrogenic_print hydrogenic_print.cpp)
target_link_libraries(hydrogenic_print mudiraclib)

add_executable(map_E_basins map_E_basins.cpp)
target_link_libraries(map_E_basins mudiraclib)