cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)

project(MPM2DProject)

#file(GLOB MPM_SRC
#  "src/*.h"
#  "src/*.cpp"
#)

add_library(mdm2dp
    include/MPMProcess.hpp
    include/MPMOutputVTK.hpp
    include/MPM_Particle.hpp
    include/MPM_GridNode.hpp
    include/MPM_TimeTracker.hpp
    include/MPM_GridElement.hpp
    include/MPM_Read.hpp
    include/MPM_Material.hpp
    src/MPMProcess.cpp
    src/MPMOutputVTK.cpp
    src/MPM_Particle.cpp
    src/MPM_GridNode.cpp
    src/MPM_GridElement.cpp
    ace_mate/SmallStrainHookePlaneStrain2D.cpp
    ace_mate/SmallStrainHookePlaneStress2D.cpp
    ace_mate/AceMaterials.hpp
)


target_include_directories(
  mdm2dp
  PUBLIC
    "include"
    "ace_mate"
)

add_executable(MPM2D
  main.cpp
)

target_link_libraries(MPM2D mdm2dp)
