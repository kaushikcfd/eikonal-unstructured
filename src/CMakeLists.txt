cmake_minimum_required(VERSION 3.4)

project(gmsh_interpreter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp Mesh2D.cpp Mesh2D.h Node.cpp Node.h Element.cpp Element.h)

add_executable(gmsh_interpreter ${SOURCE_FILES})
