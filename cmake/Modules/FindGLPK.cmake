find_path(GLPK_INCLUDE_DIR glpk.h
    PATHS
    C:/ProgramData/chocolatey/lib/glpk/tools/glpk-5.0/src/
)
find_library(GLPK_LIBRARY NAMES glpk
    PATHS
    C:/ProgramData/chocolatey/lib/glpk/tools/glpk-5.0/w64/
)

# Handle finding status with CMake standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK DEFAULT_MSG GLPK_LIBRARY GLPK_INCLUDE_DIR)

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)