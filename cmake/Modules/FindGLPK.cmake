find_path(GLPK_INCLUDE_DIR glpk.h
    PATHS
    /mingw64/include         # Windows MSYS2 path for headers
    D:/a/_temp/msys64/mingw64/include
)
find_library(GLPK_LIBRARY NAMES glpk
    PATHS
    /mingw64/lib             # Windows MSYS2 path for libraries
    D:/a/_temp/msys64/mingw64/lib
)

if(GLPK_INCLUDE_DIR AND GLPK_LIBRARY)
    set(GLPK_FOUND TRUE)
    set(GLPK_LIBRARIES ${GLPK_LIBRARY})
    set(GLPK_INCLUDE_DIRS ${GLPK_INCLUDE_DIR})
else()
    set(GLPK_FOUND FALSE)
endif()

message(STATUS "GLPK_INCLUDE_DIR: ${GLPK_INCLUDE_DIR}")
message(STATUS "GLPK_LIBRARY: ${GLPK_LIBRARY}")

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)