cmake_minimum_required(VERSION 3.20)
project(Levelset2dPolygon VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Wpedantic)

# Sibling checkout dependency (relative path). Also brings in Eigen3 (found
# and linked PUBLIC by common_geometry's own CMakeLists.txt) transitively.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../common_geometry
                  ${CMAKE_CURRENT_BINARY_DIR}/common_geometry)

add_library(levelset2d_polygon INTERFACE)

target_include_directories(levelset2d_polygon
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(levelset2d_polygon INTERFACE cxx_std_17)

# INTERFACE (not PRIVATE): marching_squares.hpp exposes ns_cg::Grid2d/Polygon2d
# in public signatures, so consumers of levelset2d_polygon need
# common_geometry's (and transitively Eigen's) include paths too.
target_link_libraries(levelset2d_polygon INTERFACE common_geometry)

enable_testing()
add_subdirectory(tests)

option(LEVELSET2D_POLYGON_BUILD_EXAMPLES "Build example programs" ON)
if(LEVELSET2D_POLYGON_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
