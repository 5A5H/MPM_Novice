#ifndef _ACE_MATERIALS_HPP_
#define _ACE_MATERIALS_HPP_

#include "sms.h"


void SmallStrainHookePlaneStrain2D(double MaterialData[2],double Eps[3],double Exp[16]);

void SmallStrainHookePlaneStress2D(double MaterialData[2],double Eps[3],double Exp[16]);


#endif
