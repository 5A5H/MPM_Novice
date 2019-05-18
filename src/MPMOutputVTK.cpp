#include <MPMOutputVTK.hpp>
#include <iostream>
#include <fstream>

// Constructor of the MPMProcess
MPMOutputVTK::MPMOutputVTK()
{
  std::cout << "VTK Started" << std::endl;
}

// Destructor of the MPMProcess
MPMOutputVTK::~MPMOutputVTK()
{
  std::cout << "VTK Terminated" << std::endl;
}
