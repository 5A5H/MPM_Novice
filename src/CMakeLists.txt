target_sources(MPM2DLib
    PRIVATE
        MPM_GridElement.cpp
        MPM_GridNode.cpp
        MPM_Particle.cpp
        MPM_OutputVTK.cpp
        MPM_Process.cpp
)

target_include_directories(MPM2DLib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
