cmake_minimum_required(VERSION 2.8.3)
project(grad_spline)

add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/arc_utilities)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/sdf_tools)

project(grad_spline)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")


find_package(catkin REQUIRED COMPONENTS
roscpp
std_msgs
)

find_package(Eigen3 REQUIRED
NLopt)

catkin_package()
include_directories(
# include
  ${catkin_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/third_party/sdf_tools/include
  ${PROJECT_SOURCE_DIR}/third_party/arc_utilities/include
  ${EIGEN3_INCLUDE_DIR}
)

add_executable(bspline src/bspline.cpp)
target_link_libraries(bspline ${catkin_LIBRARIES})

# add_executable(random src/example_random.cpp src/grad_band_optimizer.cpp)
# target_link_libraries(random ${catkin_LIBRARIES} ${PROJECT_SOURCE_DIR}/third_party/sdf_tools/lib/libsdf_tools.so
# ${NLopt_LIBRARIES} /usr/local/lib/libnlopt.so)

add_executable(astar src/test_astar.cpp src/hybrid_astar.cpp)
target_link_libraries(astar ${catkin_LIBRARIES} ${PROJECT_SOURCE_DIR}/third_party/sdf_tools/lib/libsdf_tools.so)

# add_executable(state_transit src/state_transition.cpp)
# target_link_libraries(state_transit ${catkin_LIBRARIES})
# add_executable(least_square test/least_square.cpp)
# target_link_libraries(least_square ${catkin_LIBRARIES})
