add_executable(kappa_a kappa_a.cpp)
target_link_libraries(kappa_a mudiraclib)

add_executable(calc_state calc_state.cpp)
target_link_libraries(calc_state mudiraclib)

add_executable(size_effect size_effect.cpp)
target_link_libraries(size_effect mudiraclib)
