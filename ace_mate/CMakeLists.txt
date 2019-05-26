target_sources(MPM2DLib
    PRIVATE
        MPM_AceMaterials.hpp
        sms.h
        SMSUtility.c
        SmallStrainHookePlaneStrain2D.cpp
        SmallStrainHookePlaneStress2D.cpp
)

target_include_directories(MPM2DLib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
