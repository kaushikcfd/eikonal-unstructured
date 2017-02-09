cmake_minimum_required(VERSION 3.2)

project(gmsh_interpreter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wextra -Wno-unused -pg")

set(SOURCE_FILES main.cpp main.h EikonalSolver.h EikonalSolver.cpp Mesh2D.cpp Mesh2D.h Node.cpp Node.h Element.cpp Element.h)

add_executable(gmsh_interpreter ${SOURCE_FILES})
