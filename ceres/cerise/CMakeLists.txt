CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(cerise)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(CGNUPLOT $ENV{HOME}/sources/gtl/vsv_stack/cgnuplot)

FIND_PACKAGE(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${CERES_INCLUDES} ${CGNUPLOT}/include)
LINK_DIRECTORIES(${CGNUPLOT}/lib)

add_definitions(-Wall -g -ggdb)

# helloworld
ADD_EXECUTABLE(circle_fit src/circle_fit.cpp)
TARGET_LINK_LIBRARIES(circle_fit ${CERES_LIBRARIES})

# rotations only
ADD_EXECUTABLE(optimise_rotations src/optimise_rotations.cpp
    src/rotation_errors.cpp)
TARGET_LINK_LIBRARIES(optimise_rotations ${CERES_LIBRARIES} -lcgnuplot)
