cmake_minimum_required(VERSION 3.0)
project(MuDirac)

# Define the core files as static libraries
add_subdirectory(src lib)

# Now, test executables
add_subdirectory(test)