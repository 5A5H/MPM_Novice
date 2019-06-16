// 2D Material Point Method

#include <MPM_Process.hpp>
#include <MPM_OutputVTK.hpp>
#include <MPM_Particle.hpp>
#include <MPM_GridNode.hpp>
#include <MPM_GridNodeBC.hpp>
#include <MPM_TimeTracker.hpp>
#include <MPM_GridElement.hpp>
#include <MPM_SHPQ4.hpp>
#include <MPM_Read.hpp>
#include <MPM_Material.hpp>
#include <MPM_AceMaterials.hpp>
#include <MPM_HF.hpp>

#include <ELSE_System.hpp>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h> // mkdir

//----------------------------------- Global Variables ----------------------------------------------------------------
std::vector<MPMParticle> Particle;
std::vector<MPMGridNode> GridNode;
std::vector<MPMGridElement> GridElement;
std::vector<MPMMaterial> Material;


//------------------------------------------ MAIN ---------------------------------------------------------------------
int main(int argc, char** argv)
{
    std::cout << "_____________________Welcome to ELSE!____________________" << std::endl;


    ELSE::SetupEnvironment(argc, argv);
    ELSE::ReadInputFile(ELSE::InputFileName);


    // Genrate The Time Tracker
    MPMTimeTracker MPMTimings;

//------------------------------------------- Output declaration ------------------------------------------------------
    bool ParaviewOutput = true;
    std::string ParticleOutputFile = "/Users/sash/mpm_2d/data/out/TwoParticle_Particle";
    std::string GridOutputFile = "/Users/sash/mpm_2d/data/out/TwoParticle_Grid";
    int PostFrequency = 200;


//------------------------------------------ spatial discretization ---------------------------------------------------
    std::string InputfileParticle = "two_discs_particledata.txt";
    std::string InputfileNodes = "two_discs_node.txt";
    std::string InputfileGrid = "two_discs_element.txt";

    // A log file
    std::ofstream LogFile;
    LogFile.open("Else_Files/ThisElseLog", std::ios::out);
    LogFile << "test" << std::endl;
    LogFile.close();


    //Read and Create Objects
    ReadParticle(InputfileParticle, Particle);
    ReadGridNodes(InputfileNodes, GridNode);
    ReadGridElementsQ4(InputfileGrid, GridElement);
    std::cout << "Problem Data: " << std::endl;
    std::cout << "Number of Particles    : " << Particle.size() << std::endl;
    std::cout << "Number of Grid Nodes   : " << GridNode.size() << std::endl;
    std::cout << "Number of Grid Elements: " << GridElement.size() << std::endl;


  MPMTimings.printTimeTable();
  std::cout << "_________________________ The End ________________________\n";
  MPMTimings.SetTime("Program Finish");
  return 0;
}
