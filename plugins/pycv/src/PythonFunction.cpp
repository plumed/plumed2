/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019-2023 of Toni Giorgino, Daniele Rapetti

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
#include "PythonFunction.h"

#include "plumed/core/ActionRegister.h"
#include "plumed/core/PlumedMain.h" // cite

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

#include <string>
#include <cmath>


using namespace std;
namespace py = pybind11;


namespace PLMD {
namespace pycv {

//+PLUMEDOC FUNCTION PYTHONFUNCTION
/*
Define collective variables in the Python language.

A Python module named in the `IMPORT` keyword is first imported.  The
function passed as the `FUNCTION` keyword is called at each time
step. It is assumed to receive a numpy array of shape `(N,1)`, N being
the number of arguments passed at the `ARG` keyword. See also \ref
CUSTOM.

The function should return two values: a scalar (the CV value), and
its gradient with respect to each coordinate (an array of the same
shape as the input). Not returning the gradient will prevent biasing
from working (with warnings).

Automatic differentiation and transparent compilation (including to
GPU) can be performed via Google's [JAX
library](https://github.com/google/jax): see \ref PYTHONCV.


\par Examples

The following example mimics the one in \ref CUSTOM.

\plumedfile
dAB: DISTANCE ATOMS=10,12
dAC: DISTANCE ATOMS=10,15
diff: PYTHONFUNCTION ARG=dAB,dAC IMPORT=pythonfunction FUNCTION=diff PERIODIC=NO
METAD ARG=diff WIDTH=0.1 HEIGHT=0.5 BIASFACTOR=10 PACE=100

\endplumedfile

The file `pythonfunction.py` should contain something as follows.

@code{.py}
import jax.numpy as np

def diff_f(X):
    x, y = X
    return y-x


# Just to demonstrate how auto grad is done.
# In this specific case it is just np.array([-1., 1.])

diff_grad = grad(diff_f)


# The function actually being called
def diff(x):
    return diff_f(x), diff_grad(x)

@endcode

\par See also

Use \ref PYTHONCV if you are dealing with atom coordinates directly.

Please see \ref PYTHONCV for installation and automatic
differentiation.

See \ref CUSTOM for a non-Python equivalent.


*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(PythonFunction,"PYFUNCTION")

void PythonFunction::registerKeywords( Keywords& keys ) {
  Function::registerKeywords( keys );
  keys.use("PERIODIC");
  keys.addInputKeyword("optional","ARG","scalar","the labels of the values from which the function is calculated");
  keys.add("compulsory","IMPORT","the python file to import, containing the function");
  keys.add("compulsory","CALCULATE",PYCV_DEFAULTCALCULATE,"the function to call");
  keys.add("compulsory","INIT",PYCV_DEFAULTINIT,"the function to call during the construction method of the function");
  keys.add("hidden","COMPONENTS","if provided, the function will return multiple components, with the names given");
  keys.addOutputComponent(PYCV_COMPONENTPREFIX.data(),"COMPONENTS","Each of the components output py the Python code, prefixed by py-");
  // Why is NOPBC not listed here?
}

// Everything being copied from Custom.cpp
PythonFunction::PythonFunction(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  ActionWithPython(ao) {
  try {
    py::gil_scoped_acquire gil;
    //Loading the python module
    std::string import;
    parse("IMPORT",import);
    std::string calculateFunName;
    //setting up the calculate function
    parse("CALCULATE",calculateFunName);
    log.printf("  will import %s and call function %s\n", import.c_str(),
               calculateFunName.c_str());
    // Initialize the module and function pointers
    pyModule = py::module::import(import.c_str());
    if (!py::hasattr(pyModule,calculateFunName.c_str())) {
      error("the function " + calculateFunName + " is not present in "+ import);
    }

    pyCalculate = pyModule.attr(calculateFunName.c_str());
    std::string initFunName;
    parse("INIT",initFunName);
    py::dict initDict;
    if(py::hasattr(pyModule,initFunName.c_str())) {
      log.printf("  will use %s during the initialization\n", initFunName.c_str());
      auto initFcn = pyModule.attr(initFunName.c_str());
      if (py::isinstance<py::dict>(initFcn)) {
        initDict = initFcn;
      } else {
        initDict = initFcn(this);
      }
    } else if(initFunName!=PYCV_DEFAULTINIT) {
      //If the default INIT is not preset, is not a problem
      error("the function "+ initFunName + " is not present in "+ import);
    }
    {
      std::vector<std::string> components;
      parseVector("COMPONENTS", components);
      if (components.size()>1) {
        error("Please define multiple COMPONENTS from INIT in python.");
      }
    }

    if(initDict.contains("COMPONENTS")) {
      if(initDict.contains("Value")) {
        error("The initialize dict cannot contain both \"Value\" and \"COMPONENTS\"");
      }
      if(!py::isinstance<py::dict>(initDict["COMPONENTS"])) {
        error("COMPONENTS must be a dictionary using with the name of the components as keys");
      }
      py::dict components=initDict["COMPONENTS"];
      for(auto comp: components) {
        auto settings = py::cast<py::dict>(comp.second);
        if(components.size()==1) { //a single component
          initializeValue(dynamic_cast<::PLMD::ActionWithValue&>(*this), settings);
          valueSettings(settings,getPntrToValue());
        } else {
          auto name=std::string(PYCV_COMPONENTPREFIX)
                    +"-"+py::cast<std::string>(comp.first);
          initializeComponent(dynamic_cast<::PLMD::ActionWithValue&>(*this),
                              name,
                              settings);
          valueSettings(settings,getPntrToComponent(name));
        }
      }

    } else if(initDict.contains("Value")) {
      py::dict settingsDict=initDict["Value"];
      initializeValue(dynamic_cast<::PLMD::ActionWithValue&>(*this),settingsDict);
      valueSettings(settingsDict,getPntrToValue());
    } else {
      warning("  WARNING: by defaults components periodicity is not set and component is added without derivatives - see manual\n");
      //this will crash with an error, beacuse periodicity is not explicitly set
      addValue();
    }

    log.printf("  with function : %s\n",calculateFunName.c_str());
    log<<"  Bibliography "
       <<plumed.cite(PYTHONCV_CITATION)
       <<"\n";
  } catch (const py::error_already_set &e) {
    error(e.what());
    //vdbg(e.what());
  }
}


// calculator
void PythonFunction::calculate() {
  try {
    py::gil_scoped_acquire gil;
    // Call the function
    py::object r = pyCalculate(this);
    if(getNumberOfComponents()>1) {		// MULTIPLE NAMED COMPONENTS
      calculateMultiComponent(r);
    } else { // SINGLE COMPONENT
      readReturn(r, getPntrToValue());
    }
  } catch (const py::error_already_set &e) {
    plumed_error_nested()<<"("<<getLabel()<<") caught a python exception:\n"<<e.what();
  }
}

void PythonFunction::readReturn(const py::object &r, Value* valPtr) {
// Is there more than 1 return value?
  if (py::isinstance<py::tuple>(r)||py::isinstance<py::list>(r)) {
    // 1st return value: CV
    py::list rl=r.cast<py::list>();
    pycvComm_t value = rl[0].cast<pycvComm_t>();
    valPtr->set(value);
    if (rl.size() > 1) {
      auto nargs = getNumberOfArguments();
      if(!valPtr->hasDerivatives()) {
        error(valPtr->getName()+" was declared without derivatives, but python returned with derivatives");
      }
      // 2nd return value: gradient: numpy array
      py::array_t<pycvComm_t> grad(rl[1]);
      if(grad.ndim() != 1 || grad.shape(0) != nargs) {
        log.printf("Error: wrong shape for the gradient return argument: should be (nargs=%lu), received %ld \n",
                   (unsigned long) nargs, grad.shape(0));
        error("PYFUNCTION returned wrong gradient shape error");
      }

      // To optimize, see "direct access"
      // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
      for(size_t i=0; i<nargs; i++) {
        valPtr->setDerivative(i,grad.at(i));
      }
    } else if (valPtr->hasDerivatives()) {
      plumed_merror(valPtr->getName()+" was declared with derivatives, but python returned none");
    }

  } else {
    // Only value returned. Might be an error as well.
    log.printf(BIASING_DISABLED);
    pycvComm_t value = r.cast<pycvComm_t>();
    valPtr->set(value);
  }
}

void PythonFunction::calculateMultiComponent(py::object &r) {

  const auto nc = getNumberOfComponents();
  if (py::isinstance<py::dict>(r)) {
    py::dict dataDict = r.cast<py::dict>(); // values
    for(int i=0; i < nc; ++i) {
      auto component=getPntrToComponent(i);
      //get the without "label.prefix-"
      std::string key=component->getName().substr(
                        2 + getLabel().size()
                        +PYCV_COMPONENTPREFIX.size());
      if (dataDict.contains(key.c_str())) {
        readReturn(dataDict[key.c_str()], component);
      } else {
        error( "python did not returned " + key );
      }
    }
  } else {
    // In principle one could handle a "list" return case.
    error("Multi-components pyCVs need to return dictionaries");
  }
}

}// namespace pycv
}// namespace PLMD

