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

#include <condition_variable>
#include <cstring>
#include <dlfcn.h>
#include <exception>
#include <iostream>
#include <mutex>
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

// Once used, this library is never unloaded, as any Python extension (as the
// separate extension of the plugin was): when PLUMED runs within Python, the
// Python code can keep plumedCommunications (its types and functions, whose
// code is in this library) after PLUMED has closed this library; when pycv
// started the interpreter, libpython (loaded with this library) must stay with
// it, e.g. for the thread-specific data that it destroys when a thread that
// used it exits
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
  py::gil_scoped_acquire gil;
  auto modules = py::module_::import("sys").attr("modules");
  if(modules.contains("plumedCommunications")) {
    return;
  }
  static PyModuleDef definition{};
  auto m = py::module_::create_extension_module("plumedCommunications", nullptr, &definition);
  PLMD::pycv::definePlumedCommunications(m);
  modules["plumedCommunications"] = m;
}

// The Python interpreter of pycv, shared by all its actions in all the threads.
// It is the interpreter of the Python program that runs PLUMED, which pycv
// borrows and never finalizes, or one that pycv starts for the first lease and
// finalizes after the last one (or at exit). The interpreter started by pycv
// lives in a thread of its own (the owner), which initializes it, releases the
// GIL and, when asked, takes it back and finalizes it: the other threads hold
// the GIL only while they call Python, so that any thread can use it (e.g. the
// thread-MPI ranks, or Plumed objects in different threads), and Python is
// initialized and finalized in the same thread whichever thread uses it last.
class PythonRuntime {
public:
  PythonRuntime()=default;
  PythonRuntime(const PythonRuntime&)=delete;
  PythonRuntime& operator=(const PythonRuntime&)=delete;
  ~PythonRuntime() {
    // at exit, as with the static pybind11::scoped_interpreter used before
    std::lock_guard<std::mutex> lk(mutex);
    if(owner.joinable()) {
      stopOwner();
    }
  }
  // Returns true if the interpreter has been started for this lease
  bool acquire() {
    // new leases wait for the start and the finalization of an interpreter
    std::lock_guard<std::mutex> lk(mutex);
    // under the lock: Py_GetVersion writes a static buffer at each call
    checkPythonVersion();
    keepThisLibraryLoaded();
    bool started=false;
    if(leases==0 && !Py_IsInitialized()) {
      startOwner();
      started=true;
    }
    ++leases;
    return started;
  }
  // Returns true if the interpreter has been finalized with this lease
  bool release() {
    std::lock_guard<std::mutex> lk(mutex);
    plumed_assert(leases>0);
    --leases;
    if(leases==0 && owner.joinable()) {
      stopOwner();
      return true;
    }
    return false;
  }
private:
  enum class State {starting, running, failed, stopping};
  void startOwner() {
    state=State::starting;
    ownerError.clear();
    owner=std::thread(&PythonRuntime::runOwner,this);
    std::unique_lock<std::mutex> lk(ownerMutex);
    ownerCondition.wait(lk,[this] {return state!=State::starting;});
    if(state==State::failed) {
      lk.unlock();
      owner.join();
      plumed_error()<<"pycv cannot start the Python interpreter: "<<ownerError;
    }
  }
  void stopOwner() {
    {
      std::lock_guard<std::mutex> lk(ownerMutex);
      state=State::stopping;
    }
    ownerCondition.notify_all();
    owner.join();
  }
  void setState(State s) {
    {
      std::lock_guard<std::mutex> lk(ownerMutex);
      state=s;
    }
    ownerCondition.notify_all();
  }
  void runOwner() {
    try {
      // without the Python signal handlers: the owner never runs Python code,
      // and e.g. ctrl-c must stop the program that runs PLUMED, as without pycv
      py::initialize_interpreter(false);
    } catch(const std::exception& e) {
      ownerError=e.what();
      setState(State::failed);
      return;
    }
    PyThreadState* threadState=PyEval_SaveThread();
    setState(State::running);
    {
      std::unique_lock<std::mutex> lk(ownerMutex);
      ownerCondition.wait(lk,[this] {return state==State::stopping;});
    }
    PyEval_RestoreThread(threadState);
    py::finalize_interpreter();
  }
  // protects leases and the transitions of the owner
  std::mutex mutex;
  unsigned leases=0;
  std::thread owner;
  // between the owner and the thread that starts or stops it
  std::mutex ownerMutex;
  std::condition_variable ownerCondition;
  State state=State::starting;
  std::string ownerError;
};

PythonRuntime runtime;
}

namespace PLMD {
namespace pycv {

PlumedScopedPythonInterpreter::PlumedScopedPythonInterpreter(::PLMD::Log&outlog)
  :log(outlog) {
  if(runtime.acquire()) {
    std::cerr<< "------ initialized Python interpreter\n";
    log<< "------ initialized Python interpreter\n";
  } else {
    std::cerr << "------ Python interpreter already initializated\n";
    log << "------ Python interpreter already initializated\n";
  }
  bool failed=false;
  std::string failure;
  try {
    registerPlumedCommunications();
  } catch(const std::exception& e) {
    // the exception (maybe with Python objects) is destroyed before finalization
    failed=true;
    failure=e.what();
  }
  if(failed) {
    // the destructor will not run
    if(runtime.release()) {
      log << "------ Python interpreter finalized\n";
    }
    plumed_error()<<"pycv cannot define the module plumedCommunications: "<<failure;
  }
}

PlumedScopedPythonInterpreter::~PlumedScopedPythonInterpreter() {
  // after the destruction of the Python objects of the action
  if(runtime.release()) {
    std::cerr << "------ Python interpreter finalized\n";
    log << "------ Python interpreter finalized\n";
  }
}

void ActionWithPython::releasePythonObjects(const std::vector<py::object*>& objects) noexcept {
  py::gil_scoped_acquire gil;
  for(auto* object : objects) {
    *object=py::object();
  }
}

ActionWithPython::ReleaseOnFailure::ReleaseOnFailure(std::vector<py::object*> o)
  : objects(std::move(o)),
    exceptions(std::uncaught_exceptions()) {}

ActionWithPython::ReleaseOnFailure::~ReleaseOnFailure() {
  if(std::uncaught_exceptions()>exceptions) {
    releasePythonObjects(objects);
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
