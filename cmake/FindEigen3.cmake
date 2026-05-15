# FindEigen3.cmake
# Bridges a FetchContent-populated Eigen to the variables that
# Ceres' own find_package(Eigen3) expects.

if(TARGET Eigen3::Eigen)
  set(Eigen3_FOUND TRUE)
  set(EIGEN3_FOUND TRUE)
  get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
  set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
  set(EIGEN3_VERSION "3.4.0")
else()
  find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
    PATH_SUFFIXES eigen3 eigen)
  if(EIGEN3_INCLUDE_DIR)
    set(Eigen3_FOUND TRUE)
    set(EIGEN3_FOUND TRUE)
    set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
    if(NOT TARGET Eigen3::Eigen)
      add_library(Eigen3::Eigen INTERFACE IMPORTED)
      set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    endif()
  else()
    set(Eigen3_FOUND FALSE)
    set(EIGEN3_FOUND FALSE)
  endif()
endif()
