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

typedef float pycv_t;

namespace PLMD {
namespace PythonCV {

//+PLUMEDOC COLVAR PYTHONCV
/*
Define collective variables in the Python language.

The CV function should be defined in a Python file and named `cv`. It
is assumed to receive a (N,3) numpy array with the coordinates of the
atoms defined in the action. It should return a scalar (the CV value).


\par Examples

TBD


\par Notes

CVs are differentiated automatically and compiled to native code
upon PLUMED initialization.


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

  bool no_gradient=false;

  int natoms;
  bool pbc;


public:
  static void registerKeywords( Keywords& keys );
  explicit PythonCV(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(PythonCV,"PYTHONCV")

void PythonCV::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the list of atoms to be passed to the function");
  keys.add("optional","STYLE","Python types, one of NATIVE, NUMPY or JAX");
  keys.add("compulsory","IMPORT","the python file to import, containing the function");
  keys.add("optional","FUNCTION","the function to call (defaults to CV)");
  
  // Why is NOPBC not listed here?
}

PythonCV::PythonCV(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
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

  // log.printf("  some debug info here %s",py_import_file.c_str());

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
  py::list r = py_fcn(py_X);

  // 1st return value: CV
  pycv_t value = r[0].cast<pycv_t>(); 
  setValue(value);

  // 2nd return value: gradient: numpy array of (natoms, 3)
  py::array_t<pycv_t> grad(r[1]);
  
  if(grad.ndim() != 2 ||
     grad.shape(0) != natoms ||
     grad.shape(1) != 3) {
    if(!no_gradient) {
      log.printf("Error: wrong shape for the second return argument - should be (natoms,3), is %d x %d",
		  grad.shape(0), grad.shape(1));
      no_gradient=true;
    }
  }
    

  // To optimize, see "direct access" https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
  for(int i=0; i<natoms; i++) {
    if(!no_gradient) {
      Vector3d gi(grad.at(i,0),
		  grad.at(i,1),
		  grad.at(i,2));
      setAtomsDerivatives(i,gi);
    } else {
      setAtomsDerivatives(i,Vector3d(0.,0.,0.));
    }
  }

  setBoxDerivativesNoPbc();	// ??

}

}
}



