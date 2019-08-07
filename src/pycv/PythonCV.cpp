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

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <string>
#include <cmath>


using namespace std;

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

class PythonCV : public Colvar {
    
    string py_file;
    string py_function="cv";
    string py_style="NUMPY";
    int natoms;
    bool pbc;

    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;


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
  keys.add("compulsory","FILE","the python file containing the function");
  keys.add("optional","FUNCTION","the function to call (defaults to CV)");
  
  // Why is NOPBC not listed here?
}

PythonCV::PythonCV(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  natoms = atoms.size();
  
  parse("FILE",py_file);

  parse("FUNCTION",py_function);

  parse("STYLE",py_style);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;


  checkRead();

  log << "  defined via the file " << py_file << "\n";

  log<<"  Bibliography "
     <<plumed.cite(PYTHONCV_CITATION)
     <<"\n";

  addValueWithDerivatives();
  setNotPeriodic();

  requestAtoms(atoms);

  // ----------------------------------------


  Py_Initialize();
  pName = PyUnicode_DecodeFSDefault(py_file.c_str());
  /* Error checking of pName left out */
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  if (pModule != NULL) {
      pFunc = PyObject_GetAttrString(pModule, py_function.c_str());
      /* pFunc is a new reference */
      
      if (pFunc && PyCallable_Check(pFunc)) {
	  pArgs = PyTuple_New(3*natoms);
      } else {
	  if (PyErr_Occurred())
	      PyErr_Print();
	  string err("Cannot find function: ");
	  err +=  py_function;
	  error(err);
      }
  }
  else {
      PyErr_Print();
      string err("Failed to load file: ");
      err +=  py_file;
      error(err);
  }
  

}


// calculator
void PythonCV::calculate() {

    

    
  setValue(value);
  setAtomsDerivatives(0,ga);
  setAtomsDerivatives(1,gb);
  setAtomsDerivatives(2,gc);

  setBoxDerivativesNoPbc();	// ??

}

}
}



