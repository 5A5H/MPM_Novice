#include <iostream>
#include <map>
#include <string>
#include <chrono>
#include <iterator>
#include <iomanip>

class MPMTimeTracker{
  public:
    MPMTimeTracker();
    ~MPMTimeTracker();
    void SetTime(std::string TimeTag);
    void printTimeTable();
  private:
    std::vector<std::chrono::high_resolution_clock::time_point> CodeTimings;
    std::vector<std::string> CodeTimingTags;
    double DurationSeconds();
};

MPMTimeTracker::MPMTimeTracker(){}
MPMTimeTracker::~MPMTimeTracker(){}

void MPMTimeTracker::SetTime(std::string TimeTag){
  CodeTimings.push_back(std::chrono::high_resolution_clock::now());
  CodeTimingTags.push_back(TimeTag);
}

void MPMTimeTracker::printTimeTable(){
  std::cout << "---------------------Time Table:----------------------" << std::endl;
  std::cout << "|TimeTag:                |dt to last:   |dt to start:|" << std::endl;
  int iter=0;
  for(; iter < CodeTimings.size(); iter++){
  auto delta = std::chrono::duration_cast<std::chrono::microseconds>( CodeTimings[iter] - CodeTimings[iter==0? iter: iter-1] ).count();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( CodeTimings[iter] - CodeTimings[0] ).count();
  std::cout << "|  " << std::setw(20)  << std::left << CodeTimingTags[iter]  << "  |" << std::setw(10) << delta*10e-7 << " s  |" << std::setw(10) << duration*10e-7 << " s|" << std::endl;
  }
}
