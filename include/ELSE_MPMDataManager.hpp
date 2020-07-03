//ELSE Data Management Class

/*
The Data Magager Class is used to store
the pointers to all bodies by name, to all grids by name
and last but not least the materials by name
*/

#include <map>
#include <string>

using std::map;
using std::string;

namespace ELSE {
namespace MPM {

class DataManager{
  public:
    DataManager(){};
    ~DataManager(){};
  private:
    // Map for Body pointers
    //map<string ,>
    // Map for Grid pointers
    // Map for Material pointers
};

}
}
