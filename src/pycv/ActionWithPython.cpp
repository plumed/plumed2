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

#include "ActionWithPython.h"

#include "core/ActionWithValue.h"
#include "tools/DLLoader.h"

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>
#include <Python.h>

#include <cstring>
#include <dlfcn.h>
#include <iostream>
#include <thread>


namespace py = pybind11;

namespace {
auto a=PLMD::DLLoader::EnsureGlobalDLOpen(&Py_Initialize);

// This library is compiled for one Python version: it must be the one of the
// interpreter (the check of PYBIND11_MODULE, PYBIND11_CHECK_PYTHON_VERSION)
void checkPythonVersion() {
  const char* compiled=PYBIND11_TOSTRING(PY_MAJOR_VERSION) "." PYBIND11_TOSTRING(PY_MINOR_VERSION);
  const char* running=Py_GetVersion();
  const auto length=std::strlen(compiled);
  if(std::strncmp(running,compiled,length)!=0 || (running[length]>='0' && running[length]<='9')) {
    plumed_error()<<"pycv was compiled for Python "<<compiled
                  <<", it cannot be used with Python "<<running;
  }
}

// Once plumedCommunications exists, this library is never unloaded, as any Python
// extension (as the separate extension of the plugin was): when PLUMED runs within
// Python, the Python code can keep the module (its types and functions, whose code
// is in this library) after PLUMED has closed this library; when pycv started the
// interpreter, libpython (loaded with this library) must stay with it, e.g. for the
// thread-specific data that it destroys when a thread that used it exits
void keepThisLibraryLoaded() {
  Dl_info info;
  if(dladdr(reinterpret_cast<void*>(&keepThisLibraryLoaded),&info) && info.dli_fname) {
    // the handle is never closed: the library is already loaded
    dlopen(info.dli_fname,RTLD_NOW|RTLD_NOLOAD|RTLD_NODELETE);
  }
}

// The module plumedCommunications is part of this library, not a separate
// Python extension: put it in sys.modules before any user module imports it,
// both when the interpreter is started here and when PLUMED runs within Python
// (again if the interpreter has been restarted)
void registerPlumedCommunications() {
  // as in the actions: a Python host may have released the GIL
  py::gil_scoped_acquire gil;
  auto modules = py::module_::import("sys").attr("modules");
  if(modules.contains("plumedCommunications")) {
    return;
  }
  keepThisLibraryLoaded();
  static PyModuleDef definition{};
  auto m = py::module_::create_extension_module("plumedCommunications", nullptr, &definition);
  PLMD::pycv::definePlumedCommunications(m);
  modules["plumedCommunications"] = m;
}

// The thread that started the interpreter keeps the GIL until it finalizes it
std::thread::id interpreterThread;
}

namespace PLMD {
namespace pycv {


// We can only have one interpreter globally. This is less than ideal
// because CVs can interfere with each other. The whole purpose of
// this superclass is to make a singleton. Putting it in the
// constructor makes it so that the interpreter is only initialized if
// one of the PYCV actions are used.
// https://pybind11.readthedocs.io/en/stable/reference.html#_CPPv422initialize_interpreterb

unsigned PlumedScopedPythonInterpreter::use_count=0;
std::unique_ptr<py::scoped_interpreter> PlumedScopedPythonInterpreter::interpreterGuard =
  nullptr;
std::mutex PlumedScopedPythonInterpreter::interpreterMutex{};

PlumedScopedPythonInterpreter::PlumedScopedPythonInterpreter(::PLMD::Log&outlog)
  :log(outlog) {
  std::lock_guard<std::mutex> lk(interpreterMutex);
  // under the lock: Py_GetVersion writes a static buffer at each call
  checkPythonVersion();
  if(use_count==0 && Py_IsInitialized()) {
    //this addresses the "calling pycv within a python interpreter problem"
    ++use_count;
  }

  if(use_count==0) {
    std::cerr<< "------ initialized Python interpreter\n";
    log<< "------ initialized Python interpreter\n";
    interpreterGuard = std::make_unique<py::scoped_interpreter>();
    interpreterThread = std::this_thread::get_id();
  } else {
    // Another thread would wait forever for the GIL (e.g. a second thread-MPI rank)
    if(interpreterGuard && interpreterThread!=std::this_thread::get_id()) {
      plumed_error()<<"pycv cannot use its Python interpreter from more than one thread:"
                    <<" it has been started by another thread (another thread-MPI rank,"
                    <<" or a Plumed object created in another thread)";
    }
    std::cerr << "------ Python interpreter already initializated\n";
    log << "------ Python interpreter already initializated\n";
  }
  ++use_count;
  try {
    registerPlumedCommunications();
  } catch(...) {
    // the destructor will not run
    --use_count;
    if(use_count==0) {
      interpreterGuard.reset(nullptr);
    }
    throw;
  }
}

PlumedScopedPythonInterpreter::~PlumedScopedPythonInterpreter() {
  // Finalization is tricky, because it should happen AFTER the
  // destruction of ALL the python declared variables
  std::lock_guard<std::mutex> lk(interpreterMutex);
  --use_count;
  if(use_count==0) {
    interpreterGuard.reset(nullptr);
    std::cerr << "------ Python interpreter finalized\n";
    log << "------ Python interpreter finalized\n";
  }
}

ActionWithPython::ActionWithPython (const ActionOptions&ao)
  :Action(ao),guard(log) {}

void ActionWithPython::pyParseFlag(const char* key, const ::pybind11::dict &initDict, bool& returnValue) {
  parseFlag(key, returnValue);
  if(initDict.contains(key)) {
    bool defaultRet;
    keywords.getLogicalDefault(key,defaultRet);
    if (returnValue!=defaultRet) {
      error(std::string("you specified the same keyword ").append(key)+ " both in python and in the settings file");
    }
    returnValue = initDict[key].cast<bool>();
  }
}
/******************************************************************************/
//Value/components init
void initializeValue(::PLMD::ActionWithValue& action,pybind11::dict &settingsDict) {
  action.log << "  will have a single component";
  bool withDerivatives=false;
  if(settingsDict.contains("derivative")) {
    withDerivatives=settingsDict["derivative"].cast<bool>();
  }
  if(withDerivatives) {
    action.addValueWithDerivatives();
    action.log << " WITH derivatives\n";
  } else {
    action.addValue();
    action.log << " WITHOUT derivatives\n";
  }
}

void initializeComponent(::PLMD::ActionWithValue& action,const std::string&name,py::dict &settingsDict) {
  bool withDerivatives=false;
  if(settingsDict.contains("derivative")) {
    withDerivatives=settingsDict["derivative"].cast<bool>();
  }

  if(withDerivatives) {
    action.addComponentWithDerivatives(name);
    action.log << " WITH derivatives\n";
  } else {
    action.addComponent(name);
    action.log << " WITHOUT derivatives\n";
  }
}

void valueSettings(py::dict &settings, Value* valPtr) {
  if(settings.contains("period")) {
    if (settings["period"].is_none()) {
      valPtr->setNotPeriodic();
    } else {
      py::tuple t = settings["period"];
      if(t.size()!=2) {
        plumed_merror("period must have exactly 2 components");
      }
      //the ballad py::str(t[0]).cast<std::string>() is to not care about the type of input of the user
      std::string min=py::str(t[0]).cast<std::string>();
      std::string max=py::str(t[1]).cast<std::string>();
      valPtr->setDomain(min, max);
    }
  }
}

} // namespace pycv
} // namespace PLMD
