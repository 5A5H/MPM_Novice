// System Basic Functions For ELSE
#ifndef _ELSE_SYSTEM_HPP
#define _ELSE_SYSTEM_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>

#include <sys/stat.h> // mkdir()
#include <unistd.h> //rmdir()
#include <stdlib.h> // exit, EXIT_FAILURE

#include <stdlib.h>

#include <tinyxml2.h>

namespace ELSE {
  // Variables of ELSE namespace:
    // the ELSE::LogFile is defined here
    static std::ofstream LogFile;
    // the ELSE::InputFileName is defined here
    static std::string InputFileName;

  // Functions to react on errors by tinyxml
  template<typename T>
  inline void XMLErrorCheck(T *i, std::string ErrorMessage){
    if (i == nullptr) {
      std::cout << "Error: " << ErrorMessage << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  inline void XMLErrorCheck(int i, std::string ErrorMessage){
    if (i != 0) {
      std::cout << "Error: " << ErrorMessage << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Creates standrd ELSE folders if not there and opens log file
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

    // define Logfile
    LogFile.open("Else_Files/ThisElseLog", std::ios::out);
    LogFile << "_____ELSE Log File_____" << std::endl;


    return 0;
  };

  // A specialized function to read body class from input file using tinyxml
  inline int ReadBodysFromXML(tinyxml2::XMLDocument &XMLInputFile){

    LogFile << "Reading Bodies : " << std::endl;

    // first get a pointer to the root element (xml should have only one root)
    tinyxml2::XMLNode *InpRoot = XMLInputFile.FirstChild();
    XMLErrorCheck(InpRoot,"Something went wrong parsing the input file! No root Element found");

    // Read in Body definitions
    int NoInputBodies = 0;
    std::vector<std::string> InputBodyNames;
    tinyxml2::XMLElement *BodyElement = InpRoot->FirstChildElement("MPMBODY");
    XMLErrorCheck(BodyElement,"Something went wrong parsing the input file! No MPMBODY element found!");

    while (BodyElement != nullptr){

      // get name of the body
      const char * BodyNamePointer = nullptr;
      std::string BodyName;
      BodyNamePointer = BodyElement->Attribute("Name");
      if (BodyNamePointer != nullptr) BodyName = BodyNamePointer;
      XMLErrorCheck(BodyNamePointer,"Something went wrong parsing the input file! No MPMBODY Name found");

      LogFile << "    Catch Body : " << BodyName << std::endl;

      // Now read parameter of this body
      tinyxml2::XMLElement *Parameter = BodyElement->FirstChildElement("Parameter");
      XMLErrorCheck(Parameter,"Something went wrong parsing the input file! No Parameter found on Body: "+BodyName);

      while (Parameter != nullptr) {
        // get Kind of input from attribute
        const char * KindPointer = nullptr;
        std::string Kind;
        KindPointer = Parameter->Attribute("Kind");
        if (KindPointer != nullptr) Kind = KindPointer;
        XMLErrorCheck(KindPointer,"Something went wrong parsing the input file! No Kind of MPMBody "+BodyName+" found!");
        std::cout << Kind << std::endl;

        if (Kind=="CSV") {
          // then read the file path
          const char * CSVPointer = nullptr;
          std::string CSVFile;
          CSVPointer = Parameter->Attribute("File");
          if (CSVPointer != nullptr) CSVFile = CSVPointer;
          XMLErrorCheck(CSVPointer,"Something went wrong parsing the input file! No Kind of MPMBody "+BodyName+" found!");
          std::cout << CSVFile  << std::endl;

        } else {
          std::cout << "Something went wrong parsing the input file! Only Kind=\"CSV\" supported by now!" << std::endl;
        }

      // go to next parameter of this body
      Parameter = Parameter->NextSiblingElement("Parameter");
      }

      // Iterate further for other MPMBODY definitions
      BodyElement = BodyElement->NextSiblingElement("MPMBODY");
      NoInputBodies++;
    }
    std::cout << "No of Bodies from input: " << NoInputBodies << std::endl;
    if (NoInputBodies==0) {std::cout << "Something went wrong parsing the input file! No MPMBODY element found!" << std::endl; ;return -1;}
    return 0;
  }

  // A specialized function to read material class from input file using tinyxml
  inline int ReadMaterialFromXML(tinyxml2::XMLDocument &XMLInputFile){

    LogFile << "Reading Materials : " << std::endl;

    // first get a pointer to the root element (xml should have only one root)
    tinyxml2::XMLNode *InpRoot = XMLInputFile.FirstChild();
    XMLErrorCheck(InpRoot,"Something went wrong parsing the input file! No root Element found");

    // Read in Body definitions
    int NoInputMaterials = 0;
    std::vector<std::string> InputMaterialNames;
    tinyxml2::XMLElement *MaterialElement = InpRoot->FirstChildElement("MPMMaterial");
    XMLErrorCheck(MaterialElement,"Something went wrong parsing the input file! No MPMMaterial element found!");

    while (MaterialElement != nullptr){

      // get name of the body
      const char * MaterialNamePointer = nullptr;
      std::string MaterialName;
      MaterialNamePointer = MaterialElement->Attribute("Name");
      if (MaterialNamePointer != nullptr) MaterialName = MaterialNamePointer;
      XMLErrorCheck(MaterialNamePointer,"Something went wrong parsing the input file! No MPMMaterial Name found");

      LogFile << "    Catch Material : " << MaterialName << std::endl;

      // get Kind of input from attribute
      const char * KindPointer = nullptr;
      std::string Kind;
      KindPointer = MaterialElement->Attribute("Kind");
      if (KindPointer != nullptr) Kind = KindPointer;
      XMLErrorCheck(KindPointer,"Something went wrong parsing the input file! No Kind of MPMMaterial "+MaterialName+" found!");
      LogFile << "\t\t" << MaterialName << " is of kind: " << Kind << std::endl;

      // Now read parameter of this Material
      tinyxml2::XMLElement *Parameter = MaterialElement->FirstChildElement("Parameter");
      XMLErrorCheck(Parameter,"Something went wrong parsing the input file! No Parameter found on Material: "+MaterialName);

      int NoMaterialParameters = 0;
      while (Parameter != nullptr) {

        // Read in material parameter name
        const char * ParameterNamePointer = nullptr;
        std::string ParameterName;
        ParameterNamePointer = Parameter->Attribute("Name");
        if (ParameterNamePointer != nullptr) ParameterName = ParameterNamePointer;
        XMLErrorCheck(ParameterNamePointer,"Something went wrong parsing the input file! No Name of MPMMaterial "+MaterialName+" parameter found!");

        // Read in material parameter
        // 1. get type
        const char * ParameterTypePointer = nullptr;
        std::string ParameterType;
        ParameterTypePointer = Parameter->Attribute("Type");
        if (ParameterNamePointer != nullptr) ParameterType = ParameterTypePointer;
        XMLErrorCheck(ParameterTypePointer,"Something went wrong parsing the input file! No Type of MPMMaterial "+MaterialName+" parameter "+ ParameterName +"found!");

        // 2. get the value
        if (ParameterType=="double"){

          double ParameterDoubleValue;
          ParameterDoubleValue = Parameter->DoubleAttribute("Value");
          LogFile << "\t\t" << MaterialName << " Parameter " << ParameterName << " of type double: " << ParameterDoubleValue << std::endl;

        } else if (ParameterType=="int") {

          int ParameterIntegerValue;
          ParameterIntegerValue = Parameter->IntAttribute("Value");
          LogFile << "\t\t" << MaterialName << " Parameter " << ParameterName << " of type int: " << ParameterIntegerValue << std::endl;

        } else if (ParameterType=="bool"){

          bool ParameterBoolValue;
          ParameterBoolValue = Parameter->BoolAttribute("Value");
          LogFile << "\t\t" << MaterialName << " Parameter " << ParameterName << " of type bool: " << ParameterBoolValue << std::endl;

        } else {
          XMLErrorCheck(-1," Material parameter may only be of type double, int or bool");
        }


        // go to next parameter of this material
        Parameter = Parameter->NextSiblingElement("Parameter");
        NoMaterialParameters++;
      }
      LogFile << "\t\t" << MaterialName << " has : " << NoMaterialParameters << " parameters." << std::endl;


      // Iterate further for other MPMBODY definitions
      MaterialElement = MaterialElement->NextSiblingElement("MPMMaterial");
      NoInputMaterials++;
    }
    std::cout << "No of Materials from input: " << NoInputMaterials << std::endl;
    if (NoInputMaterials==0) {std::cout << "Something went wrong parsing the input file! No MPMMaterial element found!" << std::endl; ;return -1;}
    return 0;
  }

  // The main Input file processing via tinyxml
  inline int ReadInputFile(std::string &PathToInputFile){
    // convert PathToInputFile int character array as required by tinyxml2
    char InpFile[PathToInputFile.size() + 1];
    PathToInputFile.copy(InpFile, PathToInputFile.size() + 1);
    InpFile[PathToInputFile.size()] = '\0';

    tinyxml2::XMLDocument XMLInputFile;
    XMLInputFile.LoadFile( InpFile );
    LogFile << "Reading XML-Input File: " << InpFile << std::endl;
    XMLErrorCheck(XMLInputFile.ErrorID(),"Something went wrong opening the input file!");

    ReadBodysFromXML(XMLInputFile);

    ReadMaterialFromXML(XMLInputFile);



    return 0;
  };

  // A standard output for ELSE tensors
  inline void printTensor(const std::string Name,std::array<double, 6> SymmetricTensor){
    std::string EmptyString;
    for (int i=0;i<Name.length()+3;i++) EmptyString.append(" ");
    std::cout << Name << " = " << std::setw(6) << SymmetricTensor[0] << "  "<< std::setw(6) <<  SymmetricTensor[1] << "  "<< std::setw(6) <<  SymmetricTensor[2] << std::endl;
    std::cout << EmptyString   << std::setw(6) << SymmetricTensor[1] << "  "<< std::setw(6) <<  SymmetricTensor[3] << "  "<< std::setw(6) <<  SymmetricTensor[4] << std::endl;
    std::cout << EmptyString   << std::setw(6) << SymmetricTensor[2] << "  "<< std::setw(6) <<  SymmetricTensor[4] << "  "<< std::setw(6) <<  SymmetricTensor[5] << std::endl;
  };
  inline void printTensor(const std::string Name,std::array<double, 9> Tensor){
    std::string EmptyString;
    for (int i=0;i<Name.length()+3;i++) EmptyString.append(" ");
    std::cout << Name << " = " << std::setw(6) <<  Tensor[0] << "  "<< std::setw(6) <<  Tensor[1] << "  "<< std::setw(6) <<  Tensor[2] << std::endl;
    std::cout << EmptyString   << std::setw(6) <<  Tensor[3] << "  "<< std::setw(6) <<  Tensor[4] << "  "<< std::setw(6) <<  Tensor[5] << std::endl;
    std::cout << EmptyString   << std::setw(6) <<  Tensor[6] << "  "<< std::setw(6) <<  Tensor[7] << "  "<< std::setw(6) <<  Tensor[8] << std::endl;
  };
  inline void printTensor(const std::string Name,std::array<double, 3> Tensor){
    std::cout << Name << " = " << std::setw(6) <<  Tensor[0] << "  "<< std::setw(6) <<  Tensor[1] << "  "<< std::setw(6) <<  Tensor[2] << std::endl;
  };

}// end namepace else

#endif
