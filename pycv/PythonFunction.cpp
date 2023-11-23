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
#include "colvar/ActionRegister.h"
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



class PythonFunction :
  public function::Function,
  public ActionWithPython {

  string import;
  string function_name;
  size_t nargs;

  py::array_t<pycv_t, py::array::c_style> py_arg;

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
  keys.add("compulsory","FUNCTION","the function to call");

  // Why is NOPBC not listed here?
}

// Everything being copied from Custom.cpp
PythonFunction::PythonFunction(const ActionOptions&ao):
  Action(ao),
  Function(ao)
  // PLUMED_COLVAR_INIT(ao),
  // pbc(false)
{

  nargs = getNumberOfArguments();

  parse("IMPORT",import);
  parse("FUNCTION",function_name);

  addValueWithDerivatives();
  checkRead();

  log.printf("  with function : %s\n",function_name.c_str());

  log<<"  Bibliography "
     <<plumed.cite(PYTHONCV_CITATION)
     <<"\n";


  // ----------------------------------------

  // Initialize the module and function pointer
  py_module = py::module::import(import.c_str());
  py_fcn = py_module.attr(function_name.c_str());


  // ...and the coordinates array
  py_arg = py::array_t<pycv_t>(nargs);
  // ^ 2nd template argument may be py::array::c_style if needed
  // py_X_ptr = (pycv_t *) py_X.request().ptr;

}


// calculator
void PythonFunction::calculate() {

  for(size_t i=0; i<nargs; i++) {
    py_arg.mutable_at(i)=getArgument(i);
  }

  // Call the function
  py::object r = py_fcn(py_arg);

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



}
}



