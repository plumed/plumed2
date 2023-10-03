/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019-2023 of Toni Giorgino

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
#include "PythonCVInterface.h"

#include "core/PlumedMain.h"
#include "colvar/ActionRegister.h"
#include "tools/Pbc.h"

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>
#include <pybind11/operators.h>

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

#define vdbg(...) std::cerr << std::setw(4) << __LINE__ <<":" << std::setw(20)<< #__VA_ARGS__ << " " << (__VA_ARGS__) <<'\n'

namespace py = pybind11;

using std::string;
using std::vector;

namespace PLMD {
namespace pycv {

PLUMED_REGISTER_ACTION(PythonCVInterface,"PYCVINTERFACE")

void PythonCVInterface::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the list of atoms to be passed to the function");
  keys.add("compulsory","IMPORT","the python file to import, containing the function");
  keys.add("compulsory","CALCULATE","the function to call as calculate method of a CV");
  //add other callable methods
  keys.add("optional","PREPARE","the function to call as prepare() method of the CV");
  keys.add("optional","INIT","the function to call during the construction method of the CV");
  keys.add("optional","UPDATE","the function to call as update() method of the CV");
  keys.add("optional","COMPONENTS","if provided, the function will return multiple components, with the names given");
  keys.addOutputComponent("py","COMPONENTS","Each of the components output py the Python code, prefixed by py-");
  // Why is NOPBC not listed here?
}

PythonCVInterface::PythonCVInterface(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao) {
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  natoms = atoms.size();
  if(natoms==0) error("At least one atom is required");

  parse("IMPORT",import);


  parse("CALCULATE",calculate_function);
  parse("PREPARE",prepare_function);
  parse("UPDATE",update_function);
  parse("INIT",init_function);

  parseVector("COMPONENTS",components);
  ncomponents=components.size();

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  checkRead();
  log.printf("  will import %s and call function %s\n",
             import.c_str(), calculate_function.c_str());
  log.printf("  the function will receive an array of %d x 3\n",natoms);
  if(ncomponents) {
    log.printf("  it is expected to return dictionaries with %d components\n", ncomponents);
  }


  log<<"  Bibliography "
     <<plumed.cite(PYTHONCV_CITATION)
     <<"\n";

  if(ncomponents) {
    for(auto c: components) {
      auto c_pfx="py-"+c;
      addComponentWithDerivatives(c_pfx);
      componentIsNotPeriodic(c_pfx);
    }
    log<<"  WARNING: components will not have a periodicity set - see manual\n";
  } else {
    addValueWithDerivatives();
    setNotPeriodic();
  }

  requestAtoms(atoms);

  // ----------------------------------------

  // Initialize the module and function pointer

  py_module = py::module::import(import.c_str());
  py_fcn = py_module.attr(calculate_function.c_str());
  // ^ 2nd template argument may be py::array::c_style if needed
  // py_X_ptr = (pycv_t *) py_X.request().ptr;
  if (prepare_function!=PYCV_NOTIMPLEMENTED) {
    has_prepare=true;
  }
  if (update_function!=PYCV_NOTIMPLEMENTED) {
    has_update=true;
  }
  if (init_function!=PYCV_NOTIMPLEMENTED) {
    vdbg(init_function);
    auto init_fcn = py_module.attr(init_function.c_str());
    init_fcn(this);
  }
}

void PythonCVInterface::prepare() {
  if(has_prepare) {
    auto prepare_fcn = py_module.attr(prepare_function.c_str());
    py::dict prepareDict=prepare_fcn(this);
    if (prepareDict.contains("setAtomRequest")) {
      py::tuple t=prepareDict["setAtomRequest"];
      std::vector<PLMD::AtomNumber> myatoms;
      for(const auto &i:t) {
        auto at = PLMD::AtomNumber::index(i.cast<unsigned>());
        myatoms.push_back(at);
      }
      for(const auto &i:myatoms) {
        std::cout <<i.index()<<" ";
      }
      std::cout <<"\n";
      requestAtoms(myatoms);
    }
  }
}

void PythonCVInterface::update() {
  if(has_update) {
    auto update_fcn = py_module.attr(update_function.c_str());
    py::dict updateDict=update_fcn(this);
    //See what to do here
  }
}

// calculator
void PythonCVInterface::calculate() {

  if(pbc) makeWhole();

// Call the function
  py::object r = py_fcn(this);

  if(ncomponents>0) {		// MULTIPLE NAMED COMPONENTS
    calculateMultiComponent(r);
  } else {			// SINGLE COMPONENT
    calculateSingleComponent(r);
  }

}


void PythonCVInterface::calculateSingleComponent(py::object &r) {
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
    log.printf(BIASING_DISABLED);
    pycv_t value = r.cast<pycv_t>();
    setValue(value);
  }
  setBoxDerivativesNoPbc();	// ??
}


void PythonCVInterface::calculateMultiComponent(py::object &r) {
  bool dictstyle=py::isinstance<py::dict>(r);

  if(dictstyle) {
    py::dict dataDict=r.cast<py::dict>(); // values

    for(auto c: components) {

      Value *cv=getPntrToComponent("py-"+c);

      const char *cp = c.c_str();
      auto componentData = dataDict[cp].cast<py::tuple>();
      pycv_t value = componentData[0].cast<pycv_t>();
      cv->set(value);

      py::array_t<pycv_t> grad(componentData[1]);
      check_dim(grad);

      for(int i=0; i<natoms; i++) {
        Vector3d gi(grad.at(i,0),
                    grad.at(i,1),
                    grad.at(i,2));
        setAtomsDerivatives(cv,i,gi);
      }
      setBoxDerivativesNoPbc(cv);
    }
  } else {
    // In principle one could handle a "list" return case.
    error("Sorry, multi-components needs to return dictionaries");
  }
}

// Assert correct gradient shape
void PythonCVInterface::check_dim(py::array_t<pycv_t> grad) {
  if(grad.ndim() != 2 ||
      grad.shape(0) != natoms ||
      grad.shape(1) != 3) {
    log.printf("Error: wrong shape for the gradient return argument: should be (natoms=%d,3), received %ld x %ld\n",
               natoms, grad.shape(0), grad.shape(1));
    error("Python CV returned wrong gradient shape error");
  }
}
} //pycv
} //PLMD
