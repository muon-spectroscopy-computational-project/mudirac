include(FetchContent)

FetchContent_Declare(
  dlib
  GIT_REPOSITORY "https://github.com/davisking/dlib.git"
  GIT_TAG v19.24
)
add_subdirectory(dlib)