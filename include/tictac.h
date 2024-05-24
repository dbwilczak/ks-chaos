#ifndef DW__TIC_TAC_H__
#define DW__TIC_TAC_H__

#include <chrono>
#include <ostream>

inline std::chrono::time_point<std::chrono::system_clock> tic(std::ostream& out, const char message[]){
  auto start = std::chrono::system_clock::now();
  std::time_t t = std::chrono::system_clock::to_time_t(start);
  out << "\nStart " << message << ": " << std::ctime(&t) << std::endl;
  return start;
}

inline void tac(std::ostream& out, std::chrono::time_point<std::chrono::system_clock>& start){
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;  
  auto t = std::chrono::system_clock::to_time_t(end);
  out << "Computation completed: " << std::ctime(&t)
      << "Elapsed time: " << elapsed_seconds.count()/60. << " minutes (" << elapsed_seconds.count() << " seconds)\n";    
}

#endif
