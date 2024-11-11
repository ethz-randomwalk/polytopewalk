#find_package(Catch2 REQUIRED)
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

add_executable(test_dense_walk ${TESTS_DIR}/test_dense_walk.cpp)
target_link_libraries(test_dense_walk PRIVATE utils dense sparse Catch2::Catch2WithMain ${GLPK_LIBRARY})
target_include_directories(test_dense_walk PRIVATE ${GLPK_INCLUDE_DIR})

add_executable(test_sparse_walk ${TESTS_DIR}/test_sparse_walk.cpp)
target_link_libraries(test_sparse_walk PRIVATE utils dense sparse Catch2::Catch2WithMain ${GLPK_LIBRARY})
target_include_directories(test_sparse_walk PRIVATE ${GLPK_INCLUDE_DIR})

add_executable(test_fr ${TESTS_DIR}/test_fr.cpp)
target_link_libraries(test_fr PRIVATE utils dense sparse Catch2::Catch2WithMain ${GLPK_LIBRARY})
target_include_directories(test_fr PRIVATE ${GLPK_INCLUDE_DIR})

add_executable(test_weights ${TESTS_DIR}/test_weights.cpp)
target_link_libraries(test_weights PRIVATE utils dense sparse Catch2::Catch2WithMain ${GLPK_LIBRARY})
target_include_directories(test_weights PRIVATE ${GLPK_INCLUDE_DIR})

add_executable(test_init ${TESTS_DIR}/test_init.cpp)
target_link_libraries(test_init PRIVATE utils dense sparse Catch2::Catch2WithMain ${GLPK_LIBRARY})
target_include_directories(test_init PRIVATE ${GLPK_INCLUDE_DIR})