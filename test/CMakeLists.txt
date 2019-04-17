add_library(test_main STATIC test_main.cpp)

# Tests
add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils test_main mudiraclib)

add_executable(test_integrate test_integrate.cpp)
target_link_libraries(test_integrate test_main mudiraclib)

add_executable(test_hydrogenic test_hydrogenic.cpp)
target_link_libraries(test_hydrogenic test_main mudiraclib)

add_executable(test_input test_input.cpp)
target_link_libraries(test_input test_main mudiraclib)
