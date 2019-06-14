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
