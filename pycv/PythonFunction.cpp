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

#include "core/PlumedMain.h"
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "function/Function.h"

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

extern void valueSettings( pybind11::dict &r, Value* valPtr);

class PythonFunction :
  public function::Function,
  public ActionWithPython {
  static constexpr auto PYCV_DEFAULTINIT="plumedInit";
  static constexpr auto PYCV_DEFAULTCALCULATE="plumedCalculate";
  ::pybind11::module_ pyModule {};
  ::pybind11::object pyCalculate{};
  size_t nargs;

  void check_dim(py::array_t<pycv_t> grad);
public:
  explicit PythonFunction(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(PythonFunction,"PYTHONFUNCTION")

void PythonFunction::registerKeywords( Keywords& keys ) {
  Function::registerKeywords( keys );
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory","IMPORT","the python file to import, containing the function");
  keys.add("compulsory","CALCULATE","the function to call");

  // Why is NOPBC not listed here?
}

// Everything being copied from Custom.cpp
PythonFunction::PythonFunction(const ActionOptions&ao)try:
  Action(ao),
         Function(ao),
  ActionWithPython(ao) {

  nargs = getNumberOfArguments();
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
  if(initDict.contains("Value")) {
    py::dict settingsDict=initDict["Value"];
    bool withDerivatives=false;
    if(settingsDict.contains("derivative")) {
      withDerivatives=settingsDict["derivative"].cast<bool>();
      if(withDerivatives) {
        addValueWithDerivatives();
        log << " WITH derivatives\n";
      } else {
        addValue();
        log << " WITHOUT derivatives\n";
      }
      valueSettings(settingsDict,getPntrToValue());
    } else {
      warning("  WARNING: by defaults components periodicity is not set and component is added without derivatives - see manual\n");
      //this will crash with an error, beacuse periodicity is not explicitly set
      addValue();
    }
  }

  log.printf("  with function : %s\n",calculateFunName.c_str());

  log<<"  Bibliography "
     <<plumed.cite(PYTHONCV_CITATION)
     <<"\n";

} catch (const py::error_already_set &e) {
  plumed_merror(e.what());
  //vdbg(e.what());
}


// calculator
void PythonFunction::calculate() try {
  py::array_t<pycv_t, py::array::c_style> py_arg;
  for(size_t i=0; i<nargs; i++) {
    py_arg.mutable_at(i)=getArgument(i);
  }

  // Call the function
  py::object r = pyCalculate(py_arg);

  // Is there more than 1 return value?
  if(py::isinstance<py::tuple>(r)) {
    // 1st return value: CV
    py::list rl=r.cast<py::list>();
    pycv_t value = rl[0].cast<pycv_t>();
    setValue(value);

    // 2nd return value: gradient: numpy array
    py::array_t<pycv_t> grad(rl[1]);
    check_dim(grad);

    // To optimize, see "direct access"
    // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
    for(size_t i=0; i<nargs; i++) {
      setDerivative(i,grad.at(i));
    }

  } else {
    // Only value returned. Might be an error as well.
    log.printf(BIASING_DISABLED);
    pycv_t value = r.cast<pycv_t>();
    setValue(value);
  }


} catch (const py::error_already_set &e) {
  plumed_merror(e.what());
  //vdbg(e.what());
}


// Assert correct gradient shape
void PythonFunction::check_dim(py::array_t<pycv_t> grad) {
  if(grad.ndim() != 1 ||
      grad.shape(0) != nargs) {
    log.printf("Error: wrong shape for the gradient return argument: should be (nargs=%lu), received %ld \n",
               (unsigned long) nargs, grad.shape(0));
    error("Python CV returned wrong gradient shape error");
  }
}

}// namespace pycv 
}// namespace PLMD 

