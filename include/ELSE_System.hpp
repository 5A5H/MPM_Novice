// System Basic Functions For ELSE
#ifndef _ELSE_SYSTEM_HPP
#define _ELSE_SYSTEM_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <sys/stat.h> // mkdir()
#include <unistd.h> //rmdir()
#include <stdlib.h> // exit, EXIT_FAILURE

#include <stdlib.h>

#include <tinyxml2.h>

namespace ELSE {
  template<typename T>
  void XMLErrorCheck(T *i, std::string ErrorMessage){
    if (i == nullptr) std::cout << "Error: " << ErrorMessage << std::endl;
    exit(EXIT_FAILURE);
  }
  void XMLErrorCheck(int i, std::string ErrorMessage){
    if (i != 0) std::cout << "Error: " << ErrorMessage << std::endl;
    exit(EXIT_FAILURE);
  }

  static std::string InputFileName;

  inline int SetupEnvironment(int &argc, char** &argv){

    // get input file from command line arguments
    if (argc < 2) {
      std::cout << "Error: No inpur file specified " << std::endl;
      std::exit(EXIT_FAILURE);
    }
    InputFileName = argv[1];
    std::cout << "Input File: " << InputFileName << std::endl;

    // check for existing environment
    // envoronment are two folders: 1 Else_Files Else_Post

    // try to delete the environment folders. if they are there they are gone (except when they're empty) and can be freshly created afterwards. if they doesnt exist in the first place its also fine.
    int delete_Else_Files = rmdir("Else_Files");  // 0 -> directory existed an is deleted, -1 -> directory does not exist
    int delete_Else_Post  = rmdir("Else_Post");   // 0 -> directory existed an is deleted, -1 -> directory does not exist

    // create Else_Files directory
    if (mkdir("Else_Files", 0777) == -1) {
      if (errno != 17) std::cerr << "Error :  " << strerror(errno) << std::endl;
      // if errno is 17 then the folder already exists
    }

    // create Else_Post directory
    if (mkdir("Else_Post", 0777) == -1) {
      if (errno != 17) std::cerr << "Error :  " << strerror(errno) << std::endl;
      // if errno is 17 then the folder already exists
    }

    return 0;
  };

  inline int ReadInputFile(std::string &PathToInputFile){
    // convert PathToInputFile int character array as required by tinyxml2
    char InpFile[PathToInputFile.size() + 1];
    PathToInputFile.copy(InpFile, PathToInputFile.size() + 1);
    InpFile[PathToInputFile.size()] = '\0';
std::cout << "huhu1: " << std::endl;
    tinyxml2::XMLDocument XMLInputFile;
    XMLInputFile.LoadFile( InpFile );
    //XMLErrorCheck(XMLInputFile.ErrorID(),"Something went wrong opening the input file!");
std::cout << "huhu2: " << std::endl;
    // first get a pointer to the root element (xml should have only one root)
    tinyxml2::XMLNode *InpRoot = XMLInputFile.FirstChild();
    XMLErrorCheck(InpRoot,"Something went wrong parsing the input file! No root Element found");


    // Read in Body definitions
    int NoInputBodies = 0;
    tinyxml2::XMLElement *BodyElement = InpRoot->FirstChildElement("MPMBODY");
    XMLErrorCheck(BodyElement,"Something went wrong parsing the input file! No MPMBODY element found!");

    while (BodyElement != nullptr){


      // Iterate further for other MPMBODY definitions
      BodyElement = BodyElement->NextSiblingElement("MPMBODY");
      NoInputBodies++;
    }
    std::cout << "No of Bodies from input: " << NoInputBodies << std::endl;

    //XMLElement * BodyElement = pElement->FirstChildElement("MPMBODY");

    return 0;
  };

}// end namepace else

#endif
