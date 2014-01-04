cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE RelWithDebInfo)

rosbuild_init()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(CGNUPLOT $ENV{HOME}/sources/gtl/vsv_stack/cgnuplot)
# eigen 2 or 3
find_path(EIGEN_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen3
	/usr/local/include/eigen2
	/usr/local/include/eigen
	/usr/include/eigen3
	/usr/include/eigen2
	/usr/include/eigen
	/opt/ros/diamondback/stacks/geometry/eigen/include
)

# include all libs so far
include_directories(${EIGEN_INCLUDE_DIR})


FIND_PACKAGE(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${CERES_INCLUDES} ${CGNUPLOT}/include)
LINK_DIRECTORIES(${CGNUPLOT}/lib)

add_definitions(-Wall -g -ggdb -O3)


#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
rosbuild_add_executable(es_sim src/es_sim.cpp)
#target_link_libraries(example ${PROJECT_NAME})


ADD_EXECUTABLE(optimise_rotations src/optimise_rotations.cpp
    src/rotation_errors.cpp)
TARGET_LINK_LIBRARIES(optimise_rotations ${CERES_LIBRARIES} -lcgnuplot)

ADD_EXECUTABLE(optimise_mag src/optimise_mag.cpp
    src/rotation_errors.cpp)
TARGET_LINK_LIBRARIES(optimise_mag ${CERES_LIBRARIES} -lcgnuplot)



