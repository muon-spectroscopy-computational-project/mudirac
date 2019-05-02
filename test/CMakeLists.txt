add_library(test_main STATIC test_main.cpp)

# Copying data files
set(CURRENT_DATAPATH ${CMAKE_CURRENT_BINARY_DIR})
configure_file(datapath.h.in ${CMAKE_CURRENT_SOURCE_DIR}/datapath.h @ONLY)
file(COPY data/inputtest.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)

# Tests
add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils test_main mudiraclib)

add_executable(test_integrate test_integrate.cpp)
target_link_libraries(test_integrate test_main mudiraclib)

add_executable(test_hydrogenic test_hydrogenic.cpp)
target_link_libraries(test_hydrogenic test_main mudiraclib)

add_executable(test_input test_input.cpp)
target_link_libraries(test_input test_main mudiraclib)

add_executable(test_potential test_potential.cpp)
target_link_libraries(test_potential test_main mudiraclib)

add_executable(test_transforms test_transforms.cpp)
target_link_libraries(test_transforms test_main mudiraclib)

add_executable(test_atom test_atom.cpp)
target_link_libraries(test_atom test_main mudiraclib)

