/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_Stopwatch_h
#define __PLUMED_tools_Stopwatch_h

#include <string>
#include <unordered_map>
#include <iosfwd>
#include <chrono>

namespace PLMD {

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

class Stopwatch {
/// Class to store a single stopwatch.
/// Class Stopwatch contains a collection of them
  class Watch {
  public:
    std::chrono::time_point<std::chrono::high_resolution_clock> lastStart;
    long long int total = 0;
    long long int lap = 0;
    long long int max = 0;
    long long int min = 0;
    unsigned cycles = 0;
    unsigned running = 0;
    void start();
    void stop();
    void pause();
  };
  std::unordered_map<std::string,Watch> watches;
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
void Stopwatch::start() {
  start("");
}

inline
void Stopwatch::stop() {
  stop("");
}

inline
void Stopwatch::pause() {
  pause("");
}

}


#endif
