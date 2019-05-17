#include <MPMProcess.hpp>
#include <iostream>

// Constructor of the MPMProcess
MPMProcess::MPMProcess()
{
  std::cout << "Session Started" << std::endl;
}

// Destructor of the MPMProcess
MPMProcess::~MPMProcess()
{
  std::cout << "Session Terminated" << std::endl;
}
