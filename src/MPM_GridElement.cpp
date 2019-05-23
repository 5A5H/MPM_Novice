#include <MPM_GridElement.hpp>
#include <iostream>
#include <iomanip>

// Constructor
MPMGridElement::MPMGridElement(){
}
MPMGridElement::MPMGridElement(int id, int n1, int n2, int n3, int n4){
   ID = id;
   N1 = n1;
   N2 = n2;
   N3 = n3;
   N4 = n4;
}
// Destructor s
MPMGridElement::~MPMGridElement(){
}

void MPMGridElement::Report(){
  std::cout << "------------------------------------\n";
  std::cout << "This is Element :"  << std::setw(18) << ID << std::endl;
  std::cout << "My Nodes are    :"  << N1 << ",  " << N2 << ",  " << N3 << ",  " << N4 << ",  " << std::endl;
}
