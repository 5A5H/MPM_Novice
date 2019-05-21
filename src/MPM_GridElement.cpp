#include <MPM_GridElement.hpp>
#include <iostream>
#include <iomanip>

// Constructor
MPMGridElement::MPMGridElement(){
}
// Destructor s
MPMGridElement::~MPMGridElement(){
}

void MPMGridElement::Report(){
  std::cout << "------------------------------------\n";
  std::cout << "This is Element :"  << std::setw(18) << ID << std::endl;
  std::cout << "My Nodes are    :"  << N1 << N2 << N3 << N4 << std::endl;
}
