add_library(atom STATIC atom.cpp)
add_library(boundary STATIC boundary.cpp)
add_library(hydrogenic STATIC hydrogenic.cpp)
add_library(potential STATIC potential.cpp)
add_library(integrate STATIC integrate.cpp)
add_library(utils STATIC utils.cpp)
add_library(input STATIC input.cpp)

# Define interface library
add_library(mudiraclib INTERFACE)
target_link_libraries(mudiraclib INTERFACE 
                      atom boundary hydrogenic 
                      potential integrate input utils)