/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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

#include "Exception.h"
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

Notice that as of PLUMED 2.5 it is possible to use a slightly modified
interface that allow for exception safety. In practice,
one can replace a pair of calls to Stopwatch::start() and Stopwatch::stop()
with a single call to Stopwatch::startStop(). This call will return an object
that, when goes out of scope, will stop the timer.

\notice The exception safety interace is highly recommended since it allows
to make sure that stopwatches are started and stopped consistently.

For instance the following
code
\verbatim
  {
    sw.start("A");
  // any code
    sw.stop("A");
  }
\endverbatim
can be replaced with
\verbatim
  {
    auto sww=sw.startStop("A");
  // any code

  // stopwatch is stopped when sww goes out of scope
  }
\endverbatim
Similarly, Stopwatch::startPause() can be used to replace a pair of
Stopwatch::start() and Stopwatch::pause().

The older syntax (explicitly calling `Stopwatch::start()` and `Stopwatch::pause()`) is still
allowed for backward compatibility.

Notice that the returned object is of type `Stopwatch::Handler`.
You might be willing to explicitly declare a `Stopwatch::Handler` (instead of using `auto`)
when you want to conditionally start the stopwatch. For instance:
\verbatim
  {
    Stopwatch::Handler handler;
    if(you_want_to_time_this) handler=sw.startStop();
    ... do something ...
  }
  // in case it was started, the stopwatch will stop here, at the end of the block
  // in case it was not started, nothing will happen
\endverbatim

A `Stopwatch::Handler` can not be copied but it can be moved (it behaves like a unique_ptr).
Moving it explicitly allows one to transfer it to another `Stopwatch::Handler` with a different scope.
For instance, in case you want to conditionally stop the stopwatch you might use something like this:
\verbatim
  {
    Stopwatch::Handler handler;
    if(you_want_to_time_this) handler=sw.startStop();
    ... do something ...
    if(you_want_to_stop_here) auto h2=std::move(handler);
    // the previous instruction moves handler to h2 that is then destroyed, stopping the watch
    // notice that if the stop was not started it will not stop.
    ... do something else ...
  }
  // in case it is running, the stopwatch will stop here, at the end of the block
\endverbatim

Finally, notice that in order to write the timers on an output file when the
Stopwatch is destroyed, one can store a reference to a PLMD::Log by passing it
to the Stopwatch constructor.
This will make sure timers are written also in case of a premature end.
*/

class Log;

/// Return an empty string.
/// Inline static so that it can store a static variable (for quicker access)
/// without adding a unique global symbol to a library including this header file.
inline static const std::string & StopwatchEmptyString() noexcept {
  const static std::string s;
  return s;
}

class Stopwatch {

public:
/// Forward declaration
  class Watch;
/// Auxiliary class for handling exception-safe start/pause and start/stop.
  class Handler {
    Watch* watch=nullptr;
    /// stop (true) or pause (false).
    /// might be changed to an enum if clearer.
    bool stop=false;
    /// Private constructor.
    /// This is kept private to avoid misuse. Handler objects should
    /// only be created using startPause() or startStop().
    /// stop is required to know if the destructor should stop or pause the watch.
    Handler(Watch* watch,bool stop);
    /// Allows usage of private constructor
    friend class Stopwatch;
  public:
    /// Default constructor
    Handler() = default;
    /// Default copy constructor is deleted (not copyable)
    Handler(const Handler & handler) = delete;
    /// Default copy assignment is deleted (not copyable)
    Handler & operator=(const Handler & handler) = delete;
    /// Move constructor.
    Handler(Handler && handler) noexcept;
    /// Move assignment.
    Handler & operator=(Handler && handler) noexcept;
    /// Destructor either stops or pauses the watch
    ~Handler();
  };

/// Class to store a single stopwatch.
/// Class Stopwatch contains a collection of them
  class Watch {
/// Instant in time when Watch was started last time
    std::chrono::time_point<std::chrono::high_resolution_clock> lastStart;
/// Total accumulated time, in nanoseconds
    long long int total = 0;
/// Accumulated time for this lap, in nanoseconds
    long long int lap = 0;
/// Slowest lap so far, in nanoseconds
    long long int max = 0;
/// Fastest lap so far, in nanoseconds
    long long int min = 0;
/// Total number of cycles
    unsigned cycles = 0;
/// count how many times Watch was started (+1) or stopped/paused (-1).
    unsigned running = 0;
    enum class State {started, stopped, paused};
/// keep track of state
    State state = State::stopped;
/// Allows access to internal data
    friend class Stopwatch;
  public:
/// start the watch
    Watch & start();
/// stop the watch
    Watch & stop();
/// pause the watch
    Watch & pause();
/// returns a start-stop handler
    Handler startStop();
/// returns a start-pause handler
    Handler startPause();
  };

private:

/// Pointer to a log file.
/// If set, the stopwatch is logged in its destructor.
  Log*mylog=nullptr;

/// List of watches.
/// Each watch is labeled with a string.
  std::unordered_map<std::string,Watch> watches;

/// Log over stream os.
  std::ostream& log(std::ostream& os)const;

public:
// Constructor.
  explicit Stopwatch() = default;
// Constructor.
// When destructing, stopwatch is logged.
// Make sure that log survives stopwatch. Typically, it should be declared earlier, in order
// to be destroyed later.
  explicit Stopwatch(Log&log): mylog(&log) {}
// Destructor.
  ~Stopwatch();
/// Start timer named "name"
  Stopwatch& start(const std::string&name=StopwatchEmptyString());
/// Stop timer named "name"
  Stopwatch& stop(const std::string&name=StopwatchEmptyString());
/// Pause timer named "name"
  Stopwatch& pause(const std::string&name=StopwatchEmptyString());
/// Dump all timers on an ostream
  friend std::ostream& operator<<(std::ostream&,const Stopwatch&);
/// Start with exception safety, then stop.
/// Starts the Stopwatch and returns an object that, when goes out of scope,
/// stops the watch. This allows Stopwatch to be started and stopped in
/// an exception safe manner.
  Handler startStop(const std::string&name=StopwatchEmptyString());
/// Start with exception safety, then pause.
/// Starts the Stopwatch and returns an object that, when goes out of scope,
/// pauses the watch. This allows Stopwatch to be started and paused in
/// an exception safe manner.
  Handler startPause(const std::string&name=StopwatchEmptyString());
};

inline
Stopwatch::Handler::Handler(Watch* watch,bool stop) :
  watch(watch),
  stop(stop)
{
  watch->start();
}

inline
Stopwatch::Handler::~Handler() {
  if(watch) {
    if(stop) watch->stop();
    else watch->pause();
  }
}

inline
Stopwatch& Stopwatch::start(const std::string & name) {
  watches[name].start();
  return *this;
}

inline
Stopwatch& Stopwatch::stop(const std::string & name) {
  watches[name].stop();
  return *this;
}

inline
Stopwatch& Stopwatch::pause(const std::string & name) {
  watches[name].pause();
  return *this;
}

inline
Stopwatch::Handler Stopwatch::startStop(const std::string&name) {
  return watches[name].startStop();
}

inline
Stopwatch::Handler Stopwatch::startPause(const std::string&name) {
  return watches[name].startPause();
}

inline
Stopwatch::Handler::Handler(Handler && handler) noexcept :
  watch(handler.watch),
  stop(handler.stop)
{
  handler.watch=nullptr;
}

inline
Stopwatch::Handler & Stopwatch::Handler::operator=(Handler && handler) noexcept {
  if(this!=&handler) {
    if(watch) {
      if(stop) {
        try {
          watch->stop();
        } catch(...) {
// this is to avoid problems with cppcheck, given than this method is declared as
// noexcept and stop might throw in case of an internal bug
          std::terminate();
        }
      }
      else watch->pause();
    }
    watch=handler.watch;
    stop=handler.stop;
    handler.watch=nullptr;
  }
  return *this;
}

inline
Stopwatch::Watch & Stopwatch::Watch::start() {
  state=State::started;
  running++;
  lastStart=std::chrono::high_resolution_clock::now();
  return *this;
}

inline
Stopwatch::Watch & Stopwatch::Watch::stop() {
  pause();
  state=State::stopped;
  cycles++;
  total+=lap;
  if(lap>max)max=lap;
  if(min>lap || cycles==1)min=lap;
  lap=0;
  return *this;
}

inline
Stopwatch::Watch & Stopwatch::Watch::pause() {
  state=State::paused;
// In case of an internal bug (non matching start stop within the startStop or startPause interface)
// this assertion could fail in a destructor.
// If this happens, the program should crash immediately
  plumed_assert(running>0) << "Non matching start/pause or start/stop commands in a Stopwatch";
  running--;
// notice: with exception safety the following might be converted to a plain error.
// I leave it like this for now:
  if(running!=0) return *this;
  auto t=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-lastStart);
  lap+=t.count();
  return *this;
}

inline
Stopwatch::Handler Stopwatch::Watch::startStop() {
  return Handler( this,true );
}

inline
Stopwatch::Handler Stopwatch::Watch::startPause() {
  return Handler( this,false );
}


}


#endif
