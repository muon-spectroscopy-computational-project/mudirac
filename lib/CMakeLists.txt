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

find_package(Ceres QUIET)
if(NOT Ceres_FOUND)
  message(STATUS "Ceres not found, fetching from source...")
  include(FetchContent)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

  # Prepend cmake/ so our FindGlog.cmake is used by Ceres' find_package(glog)
  list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

  # Fetch Eigen (required by Ceres, header-only). Use Populate rather than
  # MakeAvailable to avoid running Eigen's full CMake, which generates a
  # broken Eigen3Config.cmake that references a Targets file only created
  # at install time.
  find_package(Eigen3 QUIET)
  if(NOT Eigen3_FOUND)
    FetchContent_Declare(eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG        3.4.0
      GIT_SHALLOW    TRUE
      SYSTEM)
    FetchContent_Populate(eigen)
    if(NOT TARGET Eigen3::Eigen)
      add_library(Eigen3::Eigen INTERFACE IMPORTED)
      set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}"
        SYSTEM TRUE)
    endif()
    set(EIGEN3_INCLUDE_DIR "${eigen_SOURCE_DIR}" CACHE PATH "" FORCE)
    set(Eigen3_FOUND TRUE)
  endif()

  # Fetch glog (required by Ceres)
  find_package(glog QUIET)
  if(NOT glog_FOUND)
    FetchContent_Declare(glog
      GIT_REPOSITORY https://github.com/google/glog.git
      GIT_TAG        v0.6.0
      GIT_SHALLOW    TRUE
      SYSTEM)
    set(WITH_GFLAGS   OFF CACHE BOOL "" FORCE)
    set(WITH_GTEST    OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(glog)
  endif()

  # Fetch Ceres — only the dense solver is needed (2pF has 2 parameters)
  FetchContent_Declare(ceres
    GIT_REPOSITORY https://github.com/ceres-solver/ceres-solver.git
    GIT_TAG        2.2.0
    GIT_SHALLOW    TRUE
    SYSTEM)
  set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(SUITESPARSE    OFF CACHE BOOL "" FORCE)
  set(LAPACK         OFF CACHE BOOL "" FORCE)
  set(MINIGLOG       OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(ceres)
endif()

target_link_libraries(optimisation Ceres::ceres experiment)



# Define interface library
add_library(mudiraclib INTERFACE)
target_link_libraries(mudiraclib INTERFACE debugtasks config experiment output
                      atom boundary state potential
                      econfigs hydrogenic transforms
                      wavefunction integrate elements input utils optimisation)