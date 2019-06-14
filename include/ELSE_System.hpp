// System Basic Functions For ELSE
#ifndef _ELSE_SYSTEM_HPP
#define _ELSE_SYSTEM_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <sys/stat.h> // mkdir()
#include <unistd.h> //rmdir()
#include <stdlib.h>

namespace ELSE {

  std::string InputFileName;

  int SetupEnvironment(int &argc, char** &argv){

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

}// end namepace else

#endif
