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
      MPMOutputVTK(std::string SimName){
        SimulationName=SimName;
        NoParticleOutputs=0;};         // Constructor
      ~MPMOutputVTK();        // Destructor

      void TestVTUParticleExport(std::string FileName, std::vector<MPMParticle> &OutParticleContainer);
      void TestVTUGridExport(
        std::string FileName,
        std::vector<MPMGridNode> &OutNodeContainer,
        std::vector<MPMGridElement> &OutElementContainer
      );

      // Add Particle Output
      void SetOutput(std::string FileName, std::vector<MPMParticle> &OutParticleContainer, std::vector<std::string> OutputStrings);
      void WriteOutput(double &t);
      void WritePVD();

      void Report();

  private:
    std::string SimulationName;


    double CuttOff(double &Input, double CutOffValue = 1e10, double FallbackValue = 0){
            if (Input*Input > CutOffValue*CutOffValue) {
        return FallbackValue;
      } else {
        return Input;
      }
    }

    int NoParticleOutputs;
    // Storage for Particle Outputs
    // Each Particle Output (set by the SetOutput command) need to store:
    // std::string NAME with complete path
    // std::vector<std::vector<MPMParticle *>> ParticleConatiner (as pointer)
    // std::vector<std::string> Names of the output variables
    std::vector<std::string> ParticleOutputFileNames;
    std::vector<std::vector<std::string>> ParticleOutputDataNames;
    std::vector<std::vector<MPMParticle> *> ParticleOutputParticleContainers;
    std::vector<int> ParticleOutputFileWriteCounter;

    //OutputFile Container : Each wirtten file has an entry with
    //                       {Filename+Path, timestamp}
    std::vector<std::vector<std::string>> OutputFileContainer;

    // Function to write the actual VTU-File
    // The Coordinates and ParticleOutputDataNames of all Particles is written to OutputFile
    void VTUParticleExport(std::string FileName, std::vector<MPMParticle> *Particles, std::vector<std::string> ParticleOutputDataNames);

};


#endif
