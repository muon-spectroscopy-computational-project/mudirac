cmake_minimum_required(VERSION 3.0)
project(MuDirac)

# The standard is C++ 11
SET(CMAKE_CXX_STANDARD 11)

# Vendor libraries
add_subdirectory(vendor)

# Define the core files as static libraries
add_subdirectory(lib)

# Main executables
add_subdirectory(src bin)

# Now, test executables
add_subdirectory(test)