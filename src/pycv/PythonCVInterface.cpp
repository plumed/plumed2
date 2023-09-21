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
#include "PythonPlumedBase.h"

#include "core/PlumedMain.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "tools/Pbc.h"

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>
#include <pybind11/operators.h>

#include <string>
#include <cmath>
#include <iostream>


namespace py = pybind11;
using std::string;
using std::vector;


namespace PLMD {
namespace pycv {

///TODO: manual "you have to specify ATOMS=something for default atoms"
///TODO: add interface to pbc
///TODO: the topology can be assumed fixed and done on the go at each run by loading the pdb in the python code
class PythonCVInterface : public Colvar,
  public PythonPlumedBase {
  static constexpr auto PYCV_NOTIMPLEMENTED="PYCV_NOTIMPLEMENTED";
  string import;
  string calculate_function;
  string prepare_function{PYCV_NOTIMPLEMENTED};

  vector<string> components;
  int ncomponents;
  int natoms;
  bool pbc=false;
  bool has_prepare{false};

  void check_dim(py::array_t<pycv_t>);
  void calculateSingleComponent(py::object &);
  void calculateMultiComponent(py::object &);

public:
  explicit PythonCVInterface(const ActionOptions&);
// active methods:
  void calculate() override;
  void prepare() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(PythonCVInterface,"PYCVINTERFACE")

void PythonCVInterface::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the list of atoms to be passed to the function");
  keys.add("compulsory","IMPORT","the python file to import, containing the function");
  keys.add("compulsory","CALCULATE","the function to call as calculate method of a CV");
  //add other callable methods
  keys.add("optional","PREPARE","the function to call as prepare method of the CV");
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
  if(! py::isinstance<py::tuple>(r)) {        // Is there more than 1 return value?
    error("Sorry, multi-components needs to return gradients too");
  }

  // 1st return value: CV dict or array
  py::list rl=r.cast<py::list>();
  bool dictstyle=py::isinstance<py::dict>(rl[0]);

  if(dictstyle) {
    py::dict vdict=rl[0].cast<py::dict>(); // values
    py::dict gdict=rl[1].cast<py::dict>(); // gradients

    for(auto c: components) {
      Value *cv=getPntrToComponent("py-"+c);

      const char *cp = c.c_str();
      pycv_t value = vdict[cp].cast<pycv_t>();
      cv->set(value);

      py::array_t<pycv_t> grad(gdict[cp]);
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

PYBIND11_EMBEDDED_MODULE(plumedCommunications, m) {
  py::class_<PLMD::pycv::PythonCVInterface>(m, "PythonCVInterface")
  .def("getStep", [](PLMD::pycv::PythonCVInterface* self) {
    return self->getStep();
  },"Returns the current step")
  .def("getPosition",&PLMD::pycv::PythonCVInterface:: getPosition,
       "Returns the vector with the position of the requested atom \"i\"th"
       " atom requested by the action",py::arg("i"))
  .def("getPositions",[](PLMD::pycv::PythonCVInterface* self) -> py::list{
    py::list atomList;
    const auto Positions = self->getPositions();
    for(const auto &p: Positions) {
      atomList.append(p);
    }
    return atomList;
  },"Returns a list of the atomic positions of the atoms requested by the action")
  .def("getPositionsArray",[](PLMD::pycv::PythonCVInterface* self) -> py::array_t<double>{
    auto nat=self->getPositions().size();
    py::array_t<double>::ShapeContainer shape({nat,3});
    py::array_t<double> atomList(shape,
      &self->getPositions()[0][0]
    );
    return atomList;
  },"Returns a numpy.array that contaisn the atomic positions of the atoms requested by the action")
  ;

  py::class_<PLMD::Vector3d>(m, "Vector3D")
  .def(py::init<>())
  .def(py::init<double,double,double>())
  .def(py::self + py::self)//tested
  .def(py::self - py::self)//tested
  .def(py::self += py::self)
  .def(py::self -= py::self)
  .def(py::self *= float())
  .def(py::self /= float())
  .def(float() * py::self)//tested
  .def(py::self * float())//tested
  .def(-py::self)//tested
  .def("modulo",&PLMD::Vector3d::modulo,
       "Returns the module of the vector")//tested
  .def("modulo2",&PLMD::Vector3d::modulo2,
       "Returns the squared module of the vector")//tested
  .def("zero",&PLMD::Vector3d::zero,
       "Set all the vector componets to zero")
  .def("__setitem__", [](PLMD::Vector3d &self, unsigned index, double val)
  { self[index] = val; })
  .def("__getitem__", [](PLMD::Vector3d &self, unsigned index)
  { return self[index]; })
  //.def_property("x", &PLMD::Vector3d::getX, &PLMD::Vector3d::setX)
  .def("__str__", [](PLMD::Vector3d &self) {
    return "["+std::to_string(self[0])+
           ", "+std::to_string(self[1])+
           ", "+std::to_string(self[2])+
           "]";
  })
  .def("toArray",[](PLMD::Vector3d &self){
    py::array_t<double> toret{3};
    auto t= toret.mutable_unchecked<1>();
    t(0)=self[0];
    t(1)=self[1];
    t(2)=self[2];
    return toret;
  },
    "Returns a copy of the vector as a numpy[float64] array with shape \"(3,)\"");
  ;
  m.def("modulo",&PLMD::modulo<3>,
        "Returns the modulo of a Vector3D");//tested
  m.def("modulo2",&PLMD::modulo2<3>,
        "Returns the squared module of a Vector3D");//tested
  m.def("crossProduct",&PLMD::crossProduct,
        "Returns the cross product of two Vector3D");
  m.def("dotProduct",&PLMD::dotProduct<3>,
        "Returns the dot product of two Vector3D");
  m.def("delta",&PLMD::delta<3>,
        "Returns the difference of two Vector3D");
}
