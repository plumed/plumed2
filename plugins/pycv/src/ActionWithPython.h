/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019 of Toni Giorgino

The pycv module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The pycv module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_pycv_ActionWithPython_h //{
#define __PLUMED_pycv_ActionWithPython_h

#include <mutex>
#include <string>

#include "plumed/core/Action.h"

#include <pybind11/embed.h> // everything needed for embedding

namespace PLMD {
class Value;
namespace pycv {

using pycvComm_t = double;
static constexpr auto PYTHONCV_CITATION = "Giorgino, (2019). PYCV: a PLUMED 2 Module Enabling the Rapid Prototyping of Collective Variables in Python. Journal of Open Source Software, 4(42), 1773. doi:10.21105/joss.01773";
static constexpr auto BIASING_DISABLED = "PYCV: Gradient was expected as a second return value but is missing. Biasing won't work\n";
static constexpr std::string_view PYCV_COMPONENTPREFIX="py";

///This class act both as a guard for the interpreter and a a case container for the python module and functions
class PlumedScopedPythonInterpreter final {
public:
  PlumedScopedPythonInterpreter(::PLMD::Log&);
  ~PlumedScopedPythonInterpreter();
private:
  ::PLMD::Log& log;
  static unsigned use_count;
  static std::unique_ptr<::pybind11::scoped_interpreter> interpreterGuard;
  static std::mutex interpreterMutex;
};

class ActionWithPython: public virtual ::PLMD::Action {
  //the guard MUST be set up before the python objects
  // (so that it can be destroyed after them)
  PlumedScopedPythonInterpreter guard;
public:
  explicit ActionWithPython (const ActionOptions&);
  ///redefinition of parse to avoid confict between plumed.dat and python options
  template<typename T>
  void pyParse(const char* key, const ::pybind11::dict &initDict, T& returnValue);
///redefinition of parseFlag to avoid confict between plumed.dat and python options
  void pyParseFlag(const char* key, const ::pybind11::dict &initDict, bool& returnValue);
};

template<typename T>
void ActionWithPython::pyParse(
  const char* key, const ::pybind11::dict &initDict, T& returnValue) {
  T initVal(returnValue);
  parse(key,returnValue);
  //this is not robust, but with no access to Action::line we cannot use Tools::findKeyword
  if(initDict.contains(key)) {
    if (returnValue != initVal) {
      error(std::string("you specified the same keyword ").append(key)+ " both in python and in the settings file");
    }
    returnValue = initDict[key].cast<T>();
  }
}

void initializeValue(::PLMD::ActionWithValue&, pybind11::dict &);
void initializeComponent(::PLMD::ActionWithValue&,const std::string&,pybind11::dict &);
void valueSettings( pybind11::dict &r, Value* valPtr);

} // namespace pycv
} // namespace PLMD
#endif //__PLUMED_pycv_ActionWithPython_h //}
