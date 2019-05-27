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

};


#endif
