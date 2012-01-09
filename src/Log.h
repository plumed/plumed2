#ifndef __PLUMED_Log_h
#define __PLUMED_Log_h

#include <cstdio>
#include <string>
#include <sstream>

namespace PLMD{

class PlumedCommunicator;

/// Class containing the log stream.
///
/// It is similar to a FILE stream. It allows a printf() function, and
/// also to write with a << operator. Moreover, it can prefix
/// lines with the "PLUMED:" prefix, useful to grep out plumed
/// log from output
class Log
{
/// Actual FILE stream
  FILE* fp;
/// Flag to know if, when Log is destructed, file should be closed
  bool toBeClosed;
  static const int buflen=10000;
  std::ostringstream oss;
  char* buffer;
/// Communicator (to write only with the first processor)
  PlumedCommunicator& comm;
public:
/// Initialize on a given communicator
  Log(PlumedCommunicator&);
  ~Log();
/// Set a file with a specific name
  void setFile(std::string str);
/// Link to an already open FILE stream
  void set(FILE*f);
  void flush();
/// Standard printf-like function
  int printf(const char*fmt,...);
  template <class T>
  friend Log& operator<<(Log&,const T &);
};

/// Write using << syntax
template <class T>
Log& operator<<(Log&log,const T &t){
  log.oss<<t;
  log.printf("%s",log.oss.str().c_str());
  log.oss.str("");
  return log;
}


}

#endif
