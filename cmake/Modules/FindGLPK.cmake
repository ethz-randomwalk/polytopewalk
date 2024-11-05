find_path(GLPK_INCLUDE_DIR glpk.h
    PATHS
    C:/ProgramData/chocolatey/lib/glpk/tools/include    # Chocolatey install path for GLPK headers
    C:/ProgramData/chocolatey/lib/glpk/tools            # Additional path just in case
)
find_library(GLPK_LIBRARY NAMES glpk
    PATHS
    C:/ProgramData/chocolatey/lib/glpk/tools/lib        # Chocolatey install path for GLPK libraries
    C:/ProgramData/chocolatey/lib/glpk/tools
)

# Handle finding status with CMake standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK DEFAULT_MSG GLPK_LIBRARY GLPK_INCLUDE_DIR)

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)