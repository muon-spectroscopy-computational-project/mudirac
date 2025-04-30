cmake_minimum_required(VERSION 3.0)
project(MuDirac VERSION 1.0.1)

# The standard is C++ 11
set(CMAKE_CXX_STANDARD 11)

# Define the core files as static libraries
add_subdirectory(lib)

# add vendor libraries
add_subdirectory(vendor)

# Main executables
add_subdirectory(src bin)

# Now, test executables
enable_testing()
add_subdirectory(test)