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

#include <ELSE_MPMMaterial.hpp>
#include <ELSE_MPMMaterialFactory.hpp>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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


    // play with new material class
    ELSE::MPM::Material Mate1("Steel");
    Mate1.addMaterialParameter("Emod",21000.0e0);
    Mate1.addMaterialParameter("nue",0.3e0);
    Mate1.addMaterialParameter("rho",1e3);
    Mate1.addMaterialParameter("integervalue",1);
    Mate1.addMaterialParameter("boolvalue",true);
    Mate1.dumpMaterialParameter(ELSE::LogFile);

    int testint = 0;
    Mate1.getMaterialParameter("integervalue",testint);
    std::cout << "Int   : " << testint << std::endl;
    double testdbl = 1.0;
    Mate1.getMaterialParameter("Emod",testdbl);
    std::cout << "Double: " << testdbl << std::endl;
    bool testbool = false;
    Mate1.getMaterialParameter("boolvalue",testbool);
    std::cout << "Bool  : " << testbool << std::endl;

    std::array<double, 6> Sig = {0,0,0,   0,0,     0};
    std::array<double, 9> F   = {1,0,0, 0,1,0, 0,0,1};
    std::map<std::string, double> MaterialHistory;
    std::map<std::string, int>    IntegerMaterialIO;
    std::map<std::string, double> DoubleMaterialIO;
    Mate1.getCauchyStress(Sig,F,MaterialHistory,IntegerMaterialIO,DoubleMaterialIO);

    // Test the material factory
    std::string MyMaterialName;
    ELSE::MPM::Material* Mate2;
    Mate2 = ELSE::MPM::CreateMaterial("UNKNOWN","Aluminium");
    //Mate2 -> getName(MyMaterialName);
    //std::cout << "Name of the Material: " << MyMaterialName << std::endl;
    Mate2 = ELSE::MPM::CreateMaterial("LinearElasticity_A","Aluminium");
    Mate2 -> getName(MyMaterialName);
    std::cout << "Name of the Material: " << MyMaterialName << std::endl;
    Mate2 -> addMaterialParameter("Emod",7800.0e0);
    Mate2 -> addMaterialParameter("nue",0.3e0);
    Mate2 -> addMaterialParameter("rho",1e3);
    Mate2 -> dumpMaterialParameter(ELSE::LogFile);
    F = {1,0,0, 0,1,0, 0,0,1};
    ELSE::printTensor("F",F);
    Mate2 -> getCauchyStress(Sig,F,MaterialHistory,IntegerMaterialIO,DoubleMaterialIO);
    ELSE::printTensor("Sig",Sig);
    F = {1.1,0.01,0, 0,0.9,0, 0,0,0.8};
    ELSE::printTensor("F",F);
    Mate2 -> getCauchyStress(Sig,F,MaterialHistory,IntegerMaterialIO,DoubleMaterialIO);
    ELSE::printTensor("Sig",Sig);


    // Genrate The Time Tracker
    MPMTimeTracker MPMTimings;

//------------------------------------------- Output declaration ------------------------------------------------------
    bool ParaviewOutput = true;
    std::string ParticleOutputFile = "/Users/sash/mpm_2d/data/out/TwoParticle_Particle";
    std::string GridOutputFile = "/Users/sash/mpm_2d/data/out/TwoParticle_Grid";
    int PostFrequency = 200;


//------------------------------------------ spatial discretization ---------------------------------------------------
    std::string InputfileParticle = "two_disks_particledata.txt";
    std::string InputfileNodes = "two_disks_node.txt";
    std::string InputfileGrid = "two_disks_element.txt";


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
