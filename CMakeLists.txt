cmake_minimum_required(VERSION 3.2)
project(mdctMyImplementationTrial)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")

set(SOURCE_FILES mdct.c)
add_executable(mdctMyImplementationTrial ${SOURCE_FILES})

