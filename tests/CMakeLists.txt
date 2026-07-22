include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(test_marching_squares test_marching_squares.cpp)
add_executable(test_polygon_extractor test_polygon_extractor.cpp)
add_executable(test_levelset2d_polygon
  test_marching_squares.cpp
  test_polygon_extractor.cpp)

foreach(target test_marching_squares test_polygon_extractor test_levelset2d_polygon)
    target_link_libraries(${target}
        PRIVATE
            levelset2d_polygon
            GTest::gtest_main
    )
    add_test(NAME ${target} COMMAND ${target})
endforeach()
