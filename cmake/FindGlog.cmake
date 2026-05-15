# FindGlog.cmake
# Bridges a FetchContent-built glog target to the variables that
# Ceres' own find_package(glog) expects (GLOG_FOUND, GLOG_INCLUDE_DIR,
# GLOG_LIBRARY). If glog::glog is not yet available as a target, falls
# back to a standard header/library search.

if(TARGET glog::glog)
  set(GLOG_FOUND TRUE)
  get_target_property(GLOG_INCLUDE_DIR glog::glog INTERFACE_INCLUDE_DIRECTORIES)
  set(GLOG_LIBRARY glog::glog)
  set(GLOG_LIBRARIES glog::glog)
else()
  find_path(GLOG_INCLUDE_DIR glog/logging.h)
  find_library(GLOG_LIBRARY NAMES glog)
  if(GLOG_INCLUDE_DIR AND GLOG_LIBRARY)
    set(GLOG_FOUND TRUE)
    set(GLOG_LIBRARIES ${GLOG_LIBRARY})
  else()
    set(GLOG_FOUND FALSE)
  endif()
endif()
