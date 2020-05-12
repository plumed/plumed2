/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_tools_DLLoader_h
#define __PLUMED_tools_DLLoader_h

#include <stack>
#include <string>

namespace PLMD {

/// \ingroup TOOLBOX
/// Class taking care of dynamic loading.
/// It contains wrappers to the dlopen() routine.
/// It is designed so that when an object of this class goes
/// out of scope all the libraries loaded by it are unloaded. In this
/// manner, loaded libraries are automatically unloaded at the end of
/// execution. Libraries are loaded with RTDL_LOCAL option, which
/// means that they are not accessible from outside. Still, if they
/// contain self-registering classes, they will register themselves
/// to the ActionRegister object.
class DLLoader {
  std::stack<void*> handles;
  std::string lastError;
/// Deleted copy constructor
  DLLoader(const DLLoader&) = delete;
/// Deleted assignment
  DLLoader&operator=(const DLLoader&) = delete;
public:
/// Default constructor
  DLLoader();
/// Cleanup
  ~DLLoader();
/// Load a library, returning its handle
  void* load(const std::string&);
/// Returns the last error in dynamic loader
  const std::string & error();
/// Returns true if the dynamic loader is available (on some systems it may not).
  static bool installed();
};

}

#endif
