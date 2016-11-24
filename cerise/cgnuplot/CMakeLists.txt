#cmake_minimum_required(VERSION 2.4.6)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE RelWithDebInfo)

#set the default path for built executables to the "bin" directory
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
#set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

add_definitions(-g -ggdb -O0)

#common commands for building c++ executables and libraries
add_library(${PROJECT_NAME} 
	src/CGnuplot.cpp
	src/CGnuplotSurfaceGrapher.cpp
	src/CShell.cpp)
target_link_libraries(${PROJECT_NAME})

add_executable(test_cgnuplot src/test.cpp)
target_link_libraries(test_cgnuplot ${PROJECT_NAME})

add_executable(test_movie src/test_movie.cpp)
target_link_libraries(test_movie ${PROJECT_NAME})

add_executable(test_surf src/test_surf.cpp)
target_link_libraries(test_surf ${PROJECT_NAME})
