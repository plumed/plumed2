#ifndef __PLUMED_Stopwatch_h
#define __PLUMED_Stopwatch_h

#include <string>
#include <map>
#include <iosfwd>

namespace PLMD{

/**
\ingroup TOOLBOX
Class implementing stopwatch to time execution.

Each instance of this class is a container which
can keep track of several named stopwatches at
the same time. Access to the stopwatches
is obtained using start(), stop(), pause() methods,
giving as a parameter the name of the specific stopwatch.
Also an empty string can be used (un-named stopwatch).
Finally, all the times can be logged using << operator

\verbatim
#include "Stopwatch.h"

int main(){
  Stopwatch sw;
  sw.start();

  sw.start("initialization");
// do initialization ...
  sw.stop("initialization");

  for(int i=0;i<100;i++){
    sw.start("loop");
// do calculation
    sw.stop("loop");
  }

  sw.stop();
  return 0;
}

\endverbatim

Using pause a stopwatch can be put on hold until
the next start:

\verbatim
#include "Stopwatch.h"

int main(){
  Stopwatch sw;
  sw.start();

  sw.start("initialization");
// do initialization ...
  sw.stop("initialization");

  for(int i=0;i<100;i++){
    sw.start("loop");
// do calculation
    sw.pause("loop");
// here goes something that we do not want to include
    sw.start("loop");
// do calculation
    sw.stop("loop");
  }

  sw.stop();
  return 0;
}

\endverbatim

*/

class Stopwatch{
/// Class to hold the value of absolute time
  class Time{
  public:
    unsigned long sec;
/// I store nanosecond so as to allow high resolution clocks
/// (even if likely time will be measured in microseconds)
    unsigned      nsec;
    Time();
    Time operator-(const Time&)const;
    const Time & operator+=(const Time&);
    operator double()const;
    static Time get();
    void reset();
  };
/// Class to store a single stopwatch.
/// Class Stopwatch contains a collection of them
  class Watch{
  public:
    Watch();
    Time total;
    Time lastStart;
    Time lap;
    Time max;
    Time min;
    unsigned cycles;
    bool running;
    bool paused;
    void start();
    void stop();
    void pause();
  };
  std::map<std::string,Watch> watches;
  std::ostream& log(std::ostream&)const;
public:
/// Start timer named "name"
  void start(const std::string&name);
  void start();
/// Stop timer named "name"
  void stop(const std::string&name);
  void stop();
/// Pause timer named "name"
  void pause(const std::string&name);
  void pause();
/// Dump all timers on an ostream
  friend std::ostream& operator<<(std::ostream&,const Stopwatch&);
};

inline
void Stopwatch::start(){
  start("");
}

inline
void Stopwatch::stop(){
  stop("");
}

inline
void Stopwatch::pause(){
  pause("");
}

}


#endif
