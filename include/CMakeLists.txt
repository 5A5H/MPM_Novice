target_sources(MPM2DLib
    PRIVATE
        MPM_GridElement.hpp
        MPM_GridNode.hpp
        MPM_Particle.hpp
        MPM_OutputVTK.hpp
        MPM_Process.hpp
        MPM_Material.hpp
        MPM_Read.hpp
        MPM_SHPQ4.hpp
        MPM_TimeTracker.hpp
        MPM_GridNodeBC.hpp
)

target_include_directories(MPM2DLib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
