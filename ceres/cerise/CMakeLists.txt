CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(cerise)

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

# helloworld
ADD_EXECUTABLE(circle_fit src/circle_fit.cpp)
TARGET_LINK_LIBRARIES(circle_fit ${CERES_LIBRARIES})

# rotations only
ADD_EXECUTABLE(optimise_rotations src/optimise_rotations.cpp
    src/states.cpp src/rotation_errors.cpp)
TARGET_LINK_LIBRARIES(optimise_rotations ${CERES_LIBRARIES} -lcgnuplot)

# all displacement, assuming rotations have been solved already
ADD_EXECUTABLE(optimise_displacement src/optimise_displacement.cpp
    src/states.cpp)
TARGET_LINK_LIBRARIES(optimise_displacement ${CERES_LIBRARIES} -lcgnuplot)
