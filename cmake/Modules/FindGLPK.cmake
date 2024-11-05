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