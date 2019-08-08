/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "core/PlumedMain.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "tools/Pbc.h"

#include "PythonCV.h"

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

#include <string>
#include <cmath>


using namespace std;
namespace py = pybind11;


namespace PLMD {
namespace PythonCV {

//+PLUMEDOC COLVAR PYTHONCV
/*
Define collective variables in the Python language.

A Python module named in the `IMPORT` keyword is first imported.  The
function passed as the `FUNCTION` keyword is called at each time
step. It is assumed to receive a numpy array of shape `(N,3)` with the
coordinates of the `ATOMS` used in the action.

The function should return two values: a scalar (the CV value) and its
gradient with respect to each coordinate (an array of the same shape
as the input).

Automatic differentiation and transparent compilation (also to GPU)
can be performed via the [JAX
library](https://jax.readthedocs.io/en/latest/): see the example
below.


\par Examples

The following input tells plumed to print the distance between atoms 1
and 4.

\plumedfile 
cv1: PYTHONCV ATOMS=1,4 IMPORT=jaxcv FUNCTION=cv1
PRINT FILE=colvar.out ARG=*
\endplumedfile

The file `jaxcv.py` should contain something as follows.

\verbatim
import jax.numpy as np
from jax import grad, jit, vmap

def dist(x):
    d = x[0,:]-x[1,:]
    d2 = np.dot(d,d)
    return np.sqrt(d2)

grad_dist = grad(dist)

def cv1(x):
    return dist(x), grad_dist(x)

\endverbatim


\par Installation

It is best to compile the variable as a stand-alone dynamic object.
It currently does not seem to work well with Conda under OSX
(Homebrew's Python 3 is ok).

 1. Make sure you have Python 3 installed
 2. Install the JAX library: `pip3 install --user jaxlib`
 3. Go to `plumed2/src/pycv` and `make PythonCV.so` (or `.dylib` on OSX)
 4. Use `LOAD` to dynamically load the action.


*/
//+ENDPLUMEDOC


  // Unfortunately we can only have one interpreter globally. This is
  // less than ideal because CVs can interfere with each other.
  static py::scoped_interpreter guard{}; // start the interpreter and keep it alive

class PythonCV : public Colvar {
    
  string style="NUMPY";
  string import;
  string function_name="cv";

  py::module py_module;
  py::object py_fcn;

  py::array_t<pycv_t, py::array::c_style> py_X;
  pycv_t *py_X_ptr;

  int natoms;
  bool pbc;

  void check_dim(py::array_t<pycv_t> grad);



public:
  explicit PythonCV(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(PythonCV,"PYTHONCV")

void PythonCV::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the list of atoms to be passed to the function");
  keys.add("optional","STYLE","Python types, one of NATIVE, NUMPY or JAX [not implemented]");
  keys.add("compulsory","IMPORT","the python file to import, containing the function");
  keys.add("compulsory","FUNCTION","the function to call (defaults to CV)");
  
  // Why is NOPBC not listed here?
}

PythonCV::PythonCV(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  natoms = atoms.size();
  
  parse("STYLE",style);
  parse("IMPORT",import);
  parse("FUNCTION",function_name);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  checkRead();

  log.printf("  will import %s and call function %s with style %s\n",
	     import.c_str(), function_name.c_str(), style.c_str()     );
  log.printf("  the function will receive an array of %d x 3\n",natoms);

  log<<"  Bibliography "
     <<plumed.cite(PYTHONCV_CITATION)
     <<"\n";

  addValueWithDerivatives();
  setNotPeriodic();

  requestAtoms(atoms);

  // ----------------------------------------

  // Initialize the module and function pointer
  py_module = py::module::import(import.c_str());
  py_fcn = py_module.attr(function_name.c_str());


  // ...and the coordinates array
  py_X = py::array_t<pycv_t, py::array::c_style>({natoms,3}); // check if optimal layout
  // py_X_ptr = (pycv_t *) py_X.request().ptr;
 
}

  
// calculator
void PythonCV::calculate() {

  // Is there a faster way to get in bulk? We could even wrap a C++ array without copying.
  // Also, it may be faster to access the pointer rather than use "at"
  for(int i=0; i<natoms; i++) {
    Vector xi=getPosition(i);
    py_X.mutable_at(i,0) = xi[0];
    py_X.mutable_at(i,1) = xi[1];
    py_X.mutable_at(i,2) = xi[2];
  }

  // Call the function
  py::object r = py_fcn(py_X);

  if(py::isinstance<py::tuple>(r)) {
    // 1st return value: CV
    py::list rl=r.cast<py::list>();
    pycv_t value = rl[0].cast<pycv_t>(); 
    setValue(value);

    // 2nd return value: gradient: numpy array of (natoms, 3)
    py::array_t<pycv_t> grad(rl[1]);
    check_dim(grad);

    // To optimize, see "direct access"
    // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
    for(int i=0; i<natoms; i++) {
      Vector3d gi(grad.at(i,0),
		  grad.at(i,1),
		  grad.at(i,2));
      setAtomsDerivatives(i,gi);
    }
    
  } else {
    log.printf("Gradient not being returned as second return value. Biasing disabled\n");
    pycv_t value = r.cast<pycv_t>(); 
    setValue(value);
  }

  setBoxDerivativesNoPbc();	// ??

}
  

  // Assert correct gradient shape
  void PythonCV::check_dim(py::array_t<pycv_t> grad) {
    if(grad.ndim() != 2 ||
       grad.shape(0) != natoms ||
       grad.shape(1) != 3) {
      log.printf("Error: wrong shape for the second return argument - should be (natoms,3), is %d x %d\n",
		 grad.shape(0), grad.shape(1));
      error("Python output shape error");
    }
  }


  
}
}



