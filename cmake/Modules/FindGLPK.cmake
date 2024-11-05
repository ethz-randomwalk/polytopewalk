find_path(GLPK_INCLUDE_DIR glpk.h
    PATHS
    /mingw64/include         # MSYS2 path for headers
    C:/msys64/mingw64/include
    D:/a/_temp/msys64/mingw64/include
)
find_library(GLPK_LIBRARY NAMES glpk
    PATHS
    /mingw64/lib             # MSYS2 path for libraries
    C:/msys64/mingw64/lib
    D:/a/_temp/msys64/mingw64/lib
)

# Handle finding status with CMake standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK DEFAULT_MSG GLPK_LIBRARY GLPK_INCLUDE_DIR)

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)