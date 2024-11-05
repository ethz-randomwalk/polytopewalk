include(FindPackageHandleStandardArgs)

# Use Eigen3Config.cmake to do most of the work
find_package(Eigen3 CONFIG QUIET)

if(Eigen3_FOUND AND NOT EIGEN3_FOUND)
  # Special case for Eigen 3.3.4 chocolately package
  find_package_handle_standard_args(
    Eigen3
    VERSION_VAR
      Eigen3_VERSION
    REQUIRED_VARS
      Eigen3_FOUND
  )
else()
  find_package_handle_standard_args(
    Eigen3
    VERSION_VAR
      EIGEN3_VERSION_STRING
    REQUIRED_VARS
      EIGEN3_FOUND
  )
endif()
