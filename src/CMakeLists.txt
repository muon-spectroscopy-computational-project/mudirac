add_executable(kappa_a kappa_a.cpp)
target_link_libraries(kappa_a mudiraclib)

add_executable(calc_state calc_state.cpp)
target_link_libraries(calc_state mudiraclib)

add_executable(arbitrary_state arbitrary_state.cpp)
target_link_libraries(arbitrary_state mudiraclib)

add_executable(size_effect size_effect.cpp)
target_link_libraries(size_effect mudiraclib)

add_executable(hydrogenic_print hydrogenic_print.cpp)
target_link_libraries(hydrogenic_print mudiraclib)

add_executable(test_rkutta test_rkutta.cpp)
target_link_libraries(test_rkutta mudiraclib)

add_executable(map_E_basins map_E_basins.cpp)
target_link_libraries(map_E_basins mudiraclib)

add_executable(calc_Etrans calc_Etrans.cpp)
target_link_libraries(calc_Etrans mudiraclib)

add_executable(nucsize nucsize.cpp)
target_link_libraries(nucsize mudiraclib)

add_executable(print_econf print_econf.cpp)
target_link_libraries(print_econf mudiraclib)

add_executable(calc_elempot calc_elempot.cpp)
target_link_libraries(calc_elempot mudiraclib)

configure_file(mudirac.in.hpp ${CMAKE_CURRENT_SOURCE_DIR}/mudirac.hpp)
add_executable(mudirac mudirac.cpp)
target_link_libraries(mudirac mudiraclib)

add_executable(single_state single_state.cpp)
target_link_libraries(single_state mudiraclib)