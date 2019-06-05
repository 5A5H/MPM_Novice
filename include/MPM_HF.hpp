#ifndef _MPM_HELPERFUNCTIONS_HPP_
#define _MPM_HELPERFUNCTIONS_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

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


#endif
