#include <string>
#include <iostream>
#include <fstream>
#include <vector>

class MPMOutputVTK {
  public:
      MPMOutputVTK();         // Constructor
      ~MPMOutputVTK();        // Destructor
      std::string FileName;
      std::ofstream vtpfile;
      int WriteVTK();

  private:
    int WriteDataArray();
};


class MPMPoints {
  public:
  private:
    double x;
    double y;
    double z;
};

class MPMUnstructuredGrid {
  public:
  private:
    std::vector<MPMPoints> XI;
};
