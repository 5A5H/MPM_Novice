cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)

project(MPM2DProject)

add_library(mdm2dp
    include/MPM_Process.hpp
    include/MPM_OutputVTK.hpp
    include/MPM_Particle.hpp
    include/MPM_GridNode.hpp
    include/MPM_GridNodeBC.hpp
    include/MPM_Material.hpp
    include/MPM_TimeTracker.hpp
    include/MPM_GridElement.hpp
    include/MPM_Read.hpp
    include/MPM_SHPQ4.hpp
    include/MPM_HF.hpp
    ace_mate/MPM_AceMaterials.hpp
    ace_mate/sms.h
    ace_mate/SMSUtility.c
    ace_mate/SmallStrainHookePlaneStrain2D.cpp
    ace_mate/SmallStrainHookePlaneStress2D.cpp
    ace_mate/STVKPlaneStress2D.cpp
    ace_mate/STVKPlaneStrain2D.cpp
    src/MPM_Process.cpp
    src/MPM_OutputVTK.cpp
    src/MPM_Particle.cpp
    src/MPM_GridNode.cpp
    src/MPM_GridElement.cpp
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

add_executable(DiscMDM2D
  two_discs.cpp
)
target_link_libraries(DiscMDM2D mdm2dp)

add_executable(DiscMDM2Dn
  two_discs_n.cpp
)
target_link_libraries(DiscMDM2Dn mdm2dp)

add_executable(BarMDM2D
  bar.cpp
)
target_link_libraries(BarMDM2D mdm2dp)
