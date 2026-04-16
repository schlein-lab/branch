# GoogleTest via FetchContent (self-contained, no system dep).
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
    SYSTEM
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(test_delta_read test_delta_read.cpp)
target_link_libraries(test_delta_read PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_delta_read PRIVATE cxx_std_20)
gtest_discover_tests(test_delta_read)

add_executable(test_backend_vtable test_backend_vtable.cpp)
target_link_libraries(test_backend_vtable PRIVATE branch_backend GTest::gtest_main)
target_compile_features(test_backend_vtable PRIVATE cxx_std_20)
gtest_discover_tests(test_backend_vtable)

add_executable(test_features test_features.cpp)
target_link_libraries(test_features PRIVATE branch_classify GTest::gtest_main)
target_compile_features(test_features PRIVATE cxx_std_20)
gtest_discover_tests(test_features)

add_executable(test_lossless_graph test_lossless_graph.cpp)
target_link_libraries(test_lossless_graph PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_lossless_graph PRIVATE cxx_std_20)
gtest_discover_tests(test_lossless_graph)

add_executable(test_reverse_index test_reverse_index.cpp)
target_link_libraries(test_reverse_index PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_reverse_index PRIVATE cxx_std_20)
gtest_discover_tests(test_reverse_index)
