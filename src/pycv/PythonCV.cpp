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
#include "PythonPlumedBase.h"

#include "core/PlumedMain.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "tools/Pbc.h"

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

#include <string>
#include <cmath>


using namespace std;
namespace py = pybind11;


namespace PLMD {
namespace pycv {

//+PLUMEDOC COLVAR PYTHONCV
/*
Define collective variables in the Python language.

A Python module named in the `IMPORT` keyword is first imported.  The
function passed as the `FUNCTION` keyword is called at each time
step. It is assumed to receive a numpy array of shape `(N,3)` with the
coordinates of the `ATOMS` listed in the action.

The function should return two values: a scalar (the CV value), and
its gradient with respect to each coordinate (an array of the same
shape as the input). Not returning the gradient will prevent biasing
from working (with warnings).

Automatic differentiation and transparent compilation (including to
GPU) can be performed via Google's [JAX
library](https://github.com/google/jax): see the example below.


\par Examples

The following input tells PLUMED to print the distance between atoms 1
and 4.

\plumedfile
cv1: PYTHONCV ATOMS=1,4 IMPORT=distcv FUNCTION=cv
PRINT FILE=colvar.out ARG=*

\endplumedfile

The file `distcv.py` should contain something as follows.

@code{.py}
import numpy as np

# Define the distance function
def dist_f(x):
    r = x[0,:]-x[1,:]
    d2 = np.dot(r,r)
    return np.sqrt(d2)

def grad_dist(x):
    d = dist(x)
    r = x[0,:]-x[1,:]
    g = r/d
    return np.array([g,-g])

# The CV function actually called
def cv(x):
    return dist_f(x), grad_dist(x)

@endcode



\par JAX for automatic differentiation and compilation

Automatic differentiation and transparent compilation (including to
GPU) can be performed via Google's [JAX
library](https://github.com/google/jax). In a nutshell, it's sufficient
to replace `numpy` with `jax.numpy`. See the following example.


\plumedfile
cv1: PYTHONCV ATOMS=1,2,3 IMPORT=jaxcv FUNCTION=angle
PRINT FILE=colvar.out ARG=*

\endplumedfile


And, in `jaxcv.py`...

@code{.py}
# Import the JAX library
import jax.numpy as np
from jax import grad, jit, vmap

# Implementation of the angle function
def angle_f(x):
    r1 = x[0,:]-x[1,:]
    r2 = x[2,:]-x[1,:]

    costheta = np.dot(r1,r2) / np.linalg.norm(r1) / np.linalg.norm(r2)
    theta = np.arccos(costheta)
    return theta

# Use JAX to auto-gradient it
angle_grad = grad(angle)

# The CV function actually called
def angle(x):
    return angle_f(x), angle_grad(x)
@endcode


There are however
[limitations](https://github.com/google/jax#current-gotchas), the most
notable of which is that indexed assignments such as `x[i]=y` are not
allowed, and should instead be replaced by functional equivalents such
as `x=jax.ops.index_update(x, jax.ops.index[i], y)`.




\par Installation

Linking is somewhat tricky. Make sure you have Python 3 installed. It
currently does not seem to work well with Conda under OSX (Homebrew's
Python 3 is ok).  If you are feeling lucky, this may work:

\verbatim
./configure --enable-modules=+pycv --enable-python PYTHON_BIN=python3 LDFLAGS="`python3-config --ldflags`"
\endverbatim

At run time, you may need to set the `PYTHONHOME` or other
environment libraries.

Automatic differentiation examples require the JAX library: `pip3
install jaxlib`.


It may be useful to compile the variable as a stand-alone dynamic
object.  Once in the `src/pycv` directory, try `make PythonCV.so` (or
`make PythonCV.dylib` on OSX). The compilation step *should* pick
Python-specific compilation and linker flags.  Use Plumed's \ref LOAD
action to load the generated object.



\par Possible improvements (developer)

 * DONE Also enable access to other CVs instead of coordinates?
 * Pass the full atom coordinates structure directly?
 * Multicolvar in some way?
 * Pass a subset of pairwise distances instead of coordinates?
 * Box derivatives for PBC?
 * Pass the jax array directly (check if it is being copied)
 * Benchmark
 * Provide access to Plumed data, e.g.
   * Box size
   * Topology information
   * Functions for PBC wrapping/closest image


*/
//+ENDPLUMEDOC



class PythonCV : public Colvar,
  public PythonPlumedBase {

  string style="NUMPY";
  string import;
  string function_name="cv";

  py::array_t<pycv_t, py::array::c_style> py_X;
  // pycv_t *py_X_ptr;    /* For when we want to speed up */

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
  pbc(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  natoms = atoms.size();
  if(natoms==0) error("At least one atom is required");

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
  py_X = py::array_t<pycv_t>({natoms,3});
  // ^ 2nd template argument may be py::array::c_style if needed
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

  // Is there more than 1 return value?
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
    // Only value returned. Might be an error as well.
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
    log.printf("Error: wrong shape for the gradient return argument: should be (natoms=%d,3), received %d x %d\n",
               natoms, grad.shape(0), grad.shape(1));
    error("Python CV returned wrong gradient shape error");
  }
}



}
}



