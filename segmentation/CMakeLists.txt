cmake_minimum_required(VERSION 3.10)

# set the project name
project(segmentation)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -ggdb -O3 -Wall")

find_package(OpenCV REQUIRED)

# fftw
find_path(FFTW_INCLUDE_DIR fftw3.h
    /usr/local/include
    /usr/include
    /opt/local/include
)
find_library(FFTW3_LIBRARY fftw3 PATHS
    /usr/local/lib
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    /opt/local/lib
)
message("FFTW: " ${FFTW_INCLUDE_DIR})

include_directories(include ${OpenCV_INCLUDE_DIR} ${FFTW_INCLUDE_DIR})

# add the executable
add_executable(segment_labels src/segment_labels.cpp)
target_link_libraries(segment_labels m ${OpenCV_LIBRARIES}  ${FFTW3_LIBRARY})

add_executable(segment_sm src/statemachine.cpp)
target_link_libraries(segment_sm m)


