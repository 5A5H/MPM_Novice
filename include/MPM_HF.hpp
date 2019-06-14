#ifndef _MPM_HELPERFUNCTIONS_HPP_
#define _MPM_HELPERFUNCTIONS_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <sys/stat.h> // mkdir()
#include <unistd.h> //rmdir()

namespace MPM {

  // Define a class for writing data as CVS
  class File_CVS {
  public:
    File_CVS(std::string FileName){OutputFile.open(FileName, std::ios::out);};
    ~File_CVS(){OutputFile.close();};

    // This Function let you directly access the file
    std::ofstream* strm(){return &OutputFile;};

  private:
    std::ofstream OutputFile;
  };

  int SetupEnvironment(){

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

  void StatusBar(int step, double t, double tmax, double dt){
    int maxstep = t/dt;
    if(step % 20 == 0){
      std::cout << "   Time Integration Progress : [";
      for (int i = 0;i<=(t/tmax)*24;i++) std::cout << "%";
      for (int i = 0;i<=24-(t/tmax)*24;i++) std::cout << "-";
      std::cout << "]";
      std::cout << "   Progress : " << std::setprecision(3) << std::setw(4) << std::left << (t/tmax)*100 << " % \r" << std::flush;
    }
  };

}


#endif
