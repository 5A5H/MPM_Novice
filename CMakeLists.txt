cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)

project(MPM2DProject)

add_library(mdm2dp
    include/MPMProcess.hpp
    include/MPMOutputVTK.hpp
    include/MPM_Particle.hpp
    include/MPM_GridNode.hpp
    include/MPM_TimeTracker.hpp
    include/MPM_GridElement.hpp
    src/MPMProcess.cpp
    src/MPMOutputVTK.cpp
    src/MPM_Particle.cpp
    src/MPM_GridNode.cpp
    src/MPM_GridElement.cpp
)


target_include_directories(
  mdm2dp
  PUBLIC
    "include"
)

add_executable(MPM2D
  src/main.cpp
)

target_link_libraries(MPM2D mdm2dp)
