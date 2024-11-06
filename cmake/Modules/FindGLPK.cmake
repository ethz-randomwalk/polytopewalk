find_path(GLPK_INCLUDE_DIR glpk.h
    PATHS
    D:/a/polytopewalk/polytopewalk/glpk-4.65/src
    glpk-4.65/src
)
find_library(GLPK_LIBRARY NAMES glpk
    PATHS
    D:/a/polytopewalk/polytopewalk/glpk-4.65/w64
    glpk-4.65/w64
)

# Handle finding status with CMake standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK DEFAULT_MSG GLPK_LIBRARY GLPK_INCLUDE_DIR)

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)