#ifndef _MPM_OUT_PUT_VTK_HPP_
#define _MPM_OUT_PUT_VTK_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <MPM_Particle.hpp>
#include <MPM_GridNode.hpp>
#include <MPM_GridElement.hpp>

class MPMOutputVTK {
  public:
      MPMOutputVTK();         // Constructor
      ~MPMOutputVTK();        // Destructor

      void TestVTUParticleExport(std::string FileName, std::vector<MPMParticle> &OutParticleContainer);
      void TestVTUGridExport(
        std::string FileName,
        std::vector<MPMGridNode> &OutNodeContainer,
        std::vector<MPMGridElement> &OutElementContainer
      );

  private:
    double CuttOff(double &Input, double CutOffValue = 1e10, double FallbackValue = 0){
      if (Input*Input > CutOffValue*CutOffValue) {
        return FallbackValue;
      } else {
        return Input;
      }
    }

};


#endif
