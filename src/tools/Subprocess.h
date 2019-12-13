/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019 The plumed team
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
#ifndef __PLUMED_tools_Subprocess_h
#define __PLUMED_tools_Subprocess_h

#include "OFile.h"
#include "IFile.h"
#include <string>
#include <cstdio>
#include <memory>

namespace PLMD {

/// Small class to avoid including unistd.h here
class SubprocessPid;

/**
Class managing a subprocess.

The subprocess is launched and one can interact with it through a pipe.

In order not to consume resources, it might be possible to use this syntax:

\verbatim
// at construction:
Subprocess sp;
sp.stop();

// when needed
{
  auto h=sp.contStop();
  sp<<"command\n";
  sp.flush();
  sp.getline(answer);
}
// when h goes out of scope, subprocess is stopped again.
// If an exception is raised in the block, the subprocess is stopped as well.
\endverbatim

\warning
Currently `stop` and `cont` are giving problems with some MPI implementation,
In addition, notice that the stop signal is only sent to the child process and
not to the subsequently spawn processes, so it might not work as intended.
This feature is left here but is probably no a good idea to use it.
It can be enabled with `export PLUMED_ENABLE_SIGNALS=1`.

*/
class Subprocess {
  /// Process ID.
  /// We store this rather than pid_t to avoid including <unistd.h> in this header file.
  /// This remains nullptr in the child process.
  std::unique_ptr<SubprocessPid> pid;
  /// File descriptor, parent to child
  int fpc=0;
  /// File descriptor, child to parent
  int fcp=0;
  /// File pointer, parent to child
  FILE* fppc=NULL;
  /// File pointer, child to parent
  FILE* fpcp=NULL;
  /// PLUMED file object, parent to child.
  /// Used to simplify formatting
  OFile parent_to_child;
  /// PLUMED file object, child to parent.
  /// Used to simplify formatting
  IFile child_to_parent;
public:
  /// Class used to cont/stop a Subprocess in an exception safe manner.
  class Handler {
    Subprocess* sp=nullptr;
    /// Private constructor.
    /// Only to be called by Subprocess::contStop()
    explicit Handler(Subprocess* sp) noexcept;
    friend class Subprocess;
  public:
    /// Default constructor
    Handler() = default;
    /// Destructor stops the subprocess.
    ~Handler();
    /// Default copy constructor is deleted (not copyable)
    Handler(const Handler &) = delete;
    /// Default copy assignment is deleted (not copyable)
    Handler & operator=(const Handler & handler) = delete;
    /// Move constructor.
    Handler(Handler &&) noexcept;
    /// Move assignemnt.
    Handler & operator=(Handler && handler) noexcept;
  };
/// Constructor with a command line.
  explicit Subprocess(const std::string & cmd);
/// Destructor
  ~Subprocess();
/// Flush communication to process.
  void flush();
/// Check if subprocess facilities are available.
/// If it returns false, any call to Subprocess constructor will raise an exception.
  static bool available() noexcept;
/// Get a line from the subprocess.
  Subprocess & getline(std::string &);
/// Write something to the subprocess.
  template <class T> friend Subprocess& operator<<(Subprocess& ep,const T &t);
/// Send a SIGCONT to the subprocess.
/// Better used through contStop() method.
  void cont() noexcept;
/// Send a SIGSTOP to the subprocess.
/// Better used through contStop() method.
  void stop() noexcept;
/// Returns a handler to temporarily resume the process.
  Handler contStop() noexcept {
    return Handler(this);
  }
};

template <class T>
Subprocess& operator<<(Subprocess& ep,const T &t) {
  ep.parent_to_child<<t;
  return ep;
}

}

#endif

