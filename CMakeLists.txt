cmake_minimum_required(VERSION 3.15)
project(polytopewalk)

set(CMAKE_CXX_STANDARD 11)
if(UNIX)
    # For Unix-like systems, add the -fPIC option
    set(CMAKE_CXX_FLAGS "-O2 -fPIC")
elseif(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_MATH_DEFINES")
    # For Windows, specific Windows options can be set here
    get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
    foreach(dir ${dirs})
      message(STATUS "dir='${dir}'")
    endforeach()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

# Set project directories
set(POLYTOPEWALK_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCES_DIR ${POLYTOPEWALK_DIR}/src)
set(TESTS_DIR ${POLYTOPEWALK_DIR}/tests/cpp)

# Documentation configuration
if (BUILD_DOCS)
    add_subdirectory("docs")
    return() # Stop processing further targets for BUILD_DOCS
endif()

# Main project dependencies (only when not building docs)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(GLPK REQUIRED)

add_subdirectory("src")
# set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

pybind11_add_module(${PROJECT_NAME} MODULE pybind11/PybindExt.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC utils dense sparse pybind11::module Eigen3::Eigen ${GLPK_LIBRARY})
target_include_directories(${PROJECT_NAME} PRIVATE ${Eigen3_INCLUDE_DIRS} ${GLPK_INCLUDE_DIR})
install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION .)

# # Windows-specific DLL installation
# if(WIN32)
#     file(GLOB DLL_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lib/*.dll")
#     install(FILES ${DLL_FILES} DESTINATION ${PROJECT_NAME})
# endif()

# if (BUILD_TESTS)
add_subdirectory("tests")
# endif()

