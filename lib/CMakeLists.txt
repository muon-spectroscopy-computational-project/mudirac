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
add_library(experiment STATIC experiment.cpp)
add_library(debugtasks STATIC debugtasks.cpp)
add_library(optimisation STATIC optimisation.cpp)

# find and link Ceres for optimisation library
find_package(Ceres QUIET)
if (NOT Ceres_FOUND)
    message(STATUS "Ceres library not found. Fetching Ceres.")
    include(FetchContent)
    # Ceres configuration
    set(MINIGLOG ON CACHE BOOL "" FORCE) # Use bundled miniglog
    set(GFLAGS OFF CACHE BOOL "" FORCE) # Disable gflags
    set(SUITESPARSE OFF CACHE BOOL "" FORCE) # Disable SuiteSparse
    set(LAPACK OFF CACHE BOOL "" FORCE) # Disable LAPACK
    set(CXSPARSE OFF CACHE BOOL "" FORCE) # Disable CXSparse
    set(ACCELERATESPARSE OFF CACHE BOOL "" FORCE) # Disable AccelerateSparses
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE) # Disable building tests
    # Fetch Eigen
    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
    )
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE) # Disable building Eigen tests
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE) # Disable pkgconfig generation
    set(EIGEN_BUILD_DOCS OFF CACHE BOOL "" FORCE) # Disable building documentation
    FetchContent_MakeAvailable(Eigen)

    # Fetch Ceres solver
    FetchContent_Declare(
        Ceres
        GIT_REPOSITORY https://ceres-solver.googlesource.com/ceres-solver
        GIT_TAG 2.1.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Ceres)

    message(STATUS "Ceres library fetched and intalled")
else(
    message(STATUS "Found system installed ceres solver")
)
endif()

    target_link_libraries(optimisation Ceres::ceres experiment)


# Define interface library
add_library(mudiraclib INTERFACE)
target_link_libraries(mudiraclib INTERFACE debugtasks config experiment output
                      atom boundary state potential
                      econfigs hydrogenic transforms
                      wavefunction integrate elements input utils optimisation)