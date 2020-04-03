add_library(atom STATIC atom.cpp)
add_library(boundary STATIC boundary.cpp)
add_library(state STATIC state.cpp)
add_library(econfigs STATIC econfigs.cpp)
add_library(hydrogenic STATIC hydrogenic.cpp)
add_library(potential STATIC potential.cpp)
add_library(transforms STATIC transforms.cpp)
add_library(integrate STATIC integrate.cpp)
add_library(utils STATIC utils.cpp)
add_library(elements STATIC elements.cpp)
add_library(input STATIC input.cpp)
add_library(output STATIC output.cpp)
add_library(wavefunction STATIC wavefunction.cpp)
add_library(config STATIC config.cpp)
add_library(debugtasks STATIC debugtasks.cpp)

# Define interface library
add_library(mudiraclib INTERFACE)
target_link_libraries(mudiraclib INTERFACE debugtasks config output
                      atom boundary state potential
                      econfigs hydrogenic transforms 
                      wavefunction integrate elements input utils)