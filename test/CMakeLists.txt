add_library(test_main STATIC test_main.cpp)

# Copying data files
set(CURRENT_DATAPATH ${CMAKE_CURRENT_BINARY_DIR})
configure_file(datapath.in.h ${CMAKE_CURRENT_SOURCE_DIR}/datapath.h)
file(COPY data/inputtest.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
file(COPY data/kappa_alpha_exp.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)

# Tests
add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils test_main mudiraclib)
add_test(utils test_utils)

add_executable(test_elements test_elements.cpp)
target_link_libraries(test_elements test_main mudiraclib)
add_test(elements test_elements)

add_executable(test_econfigs test_econfigs.cpp)
target_link_libraries(test_econfigs test_main mudiraclib)
add_test(econfigs test_econfigs)

add_executable(test_integrate test_integrate.cpp)
target_link_libraries(test_integrate test_main mudiraclib)
add_test(integrate test_integrate)

add_executable(test_hydrogenic test_hydrogenic.cpp)
target_link_libraries(test_hydrogenic test_main mudiraclib)
add_test(hydrogenic test_hydrogenic)

add_executable(test_input test_input.cpp)
target_link_libraries(test_input test_main mudiraclib)
add_test(input test_input)

add_executable(test_potential test_potential.cpp)
target_link_libraries(test_potential test_main mudiraclib)
add_test(potential test_potential)

add_executable(test_transforms test_transforms.cpp)
target_link_libraries(test_transforms test_main mudiraclib)
add_test(transforms test_transforms)

add_executable(test_atom test_atom.cpp)
target_link_libraries(test_atom test_main mudiraclib)
add_test(atom test_atom)

add_executable(test_wavefunction test_wavefunction.cpp)
target_link_libraries(test_wavefunction test_main mudiraclib)
add_test(wavefunction test_wavefunction)

add_executable(test_lines test_lines.cpp)
target_link_libraries(test_lines test_main mudiraclib)
add_test(lines test_lines)

# Removed for now as too long and unreliable. Test_lines serves a similar purpose but more efficiently.
# add_executable(test_kappaa test_kappaa.cpp)
# target_link_libraries(test_kappaa mudiraclib)
# add_test(kappaa test_kappaa)

add_custom_target(tests)
add_dependencies(tests test_utils test_elements test_econfigs test_integrate
test_hydrogenic test_input test_potential test_transforms test_atom test_wavefunction
test_lines)