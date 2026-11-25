/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023 Daniele Rapetti

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

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl_bind.h>

#include "plumed/tools/Vector.h"
#include "plumed/tools/NeighborList.h"

#include "PythonCVInterface.h"
#include "PythonFunction.h"

namespace py=pybind11;

#define defGetter(pyfun,classname, cppfun, type, description) \
  def(pyfun, [](classname* self) -> type{ \
    return self->cppfun(); }, description)

//NB: the action methods are written two times due to the diamond inheritance

PYBIND11_MODULE(plumedCommunications, m) {
  /*******************************default submodule****************************/
  py::module_ defaults = m.def_submodule("defaults", "Submodule with the default definitions");
  defaults.attr("COMPONENT") = py::dict(py::arg("period")=py::none(),py::arg("derivative")=true);
  defaults.attr("COMPONENT_NODEV") = py::dict(py::arg("period")=py::none(),py::arg("derivative")=false);

  /*************************PLMD::pycv::PythonCVInterface**********************/
  py::class_<PLMD::pycv::PythonCVInterface>(m, "PythonCVInterface")
  .def_readwrite("data",&PLMD::pycv::PythonCVInterface::dataContainer,"Return an accessible dictionary that persist along all the simulation")
  /***************************Start of Action methods***************************/
  //using "PLMD::pycv::PythonCVInterface::getStep" instead of the lambda gives compilation errors
  .defGetter("getStep",PLMD::pycv::PythonCVInterface,getStep, long int,"Returns the current step")
  .defGetter("getTime",PLMD::pycv::PythonCVInterface,getTime,double,"Return the present time")
  .defGetter("getTimeStep",PLMD::pycv::PythonCVInterface,getTimeStep,double,"Return the timestep")
  .defGetter("isRestart",PLMD::pycv::PythonCVInterface,getRestart,bool,"Return true if we are doing a restart")
  .defGetter("isExchangeStep",PLMD::pycv::PythonCVInterface,getExchangeStep,bool,"Check if we are on an exchange step")
  .def_property_readonly("label",
  [](PLMD::pycv::PythonCVInterface* self) -> std::string {
    return self->getLabel();
  },
  "returns the label")
  .def("log",[](PLMD::pycv::PythonCVInterface* self, py::object data) {
    self->log << py::str(data).cast<std::string>();
  },
  "puts a string in the PLUMED output",py::arg("s"))
  .def("lognl",[](PLMD::pycv::PythonCVInterface* self, py::object data) {
    self->log << py::str(data).cast<std::string>()<< "\n";
  },
  "puts a string in the PLUMED output (and appends a newline)",py::arg("s"))
  /****************************End of Action methods***************************/
  .def("getPosition",
  [](PLMD::pycv::PythonCVInterface* self, int i) -> py::array_t<double> {
    py::array_t<double>::ShapeContainer shape({3});
    py::array_t<double> atom(shape,&self->getPosition(i)[0]);
    return atom;
  },
  "Returns an ndarray with the position of the \"i\"th"
  " atom requested by the action",py::arg("i"))
  .def_property_readonly("nat",
  [](PLMD::pycv::PythonCVInterface* self) -> size_t {
    return self->getPositions().size();
  },
  "return the number of atoms"
                        )

  //here I can use &PLMD::pycv::PythonCVInterface::makeWhole because is not in Action
  // that is inherithed both by withValue and Atomistic
  .def("makeWhole",&PLMD::pycv::PythonCVInterface::makeWhole,"Make atoms whole, assuming they are in the proper order")
  .def("getPositions",[](PLMD::pycv::PythonCVInterface* self) -> py::array_t<double> {
    auto nat=self->getPositions().size();
    ///nat and 3 must be the same type: no warning
    py::array_t<double>::ShapeContainer shape({nat,static_cast<decltype(nat)>(3)});
    py::array_t<double> atomList(shape,
                                 &self->getPositions()[0][0]
                                );
    return atomList;
  },
  "Returns a numpy.array that contains the atomic positions of the atoms requested by the action")
  .def("getPbc",&PLMD::pycv::PythonCVInterface::getPbc,
       "returns an interface to the current pbcs")
  .def("getNeighbourList",&PLMD::pycv::PythonCVInterface::getNL,
       "returns an interface to the current Neighborlist")

//https://pybind11.readthedocs.io/en/stable/advanced/functions.html#return-value-policies
  /*.def("getAbsoluteIndexes", &PLMD::pycv::PythonCVInterface::getAbsoluteIndexes ,
  "Get the vector of absolute indexes.",
  py::return_value_policy::reference)*/
  //using unsigned: if AtomNumber changes, also this must change
  .def("absoluteIndexes",[](PLMD::pycv::PythonCVInterface* self) -> py::array_t<unsigned> {
    auto nat=self->getPositions().size();
    py::array_t<unsigned>::ShapeContainer shape({nat});
    py::array_t<unsigned> atomIndexes(shape);
    auto retAccessor = atomIndexes.mutable_unchecked<1>();
    for (decltype(nat) i=0; i < nat; ++i) {
      //at time of writing getAbsoluteIndexes returns const std::vector<AtomNumber> &
      retAccessor(i)
      =self->getAbsoluteIndexes()[i].index();
    }
    return atomIndexes;
  },
  "Get the vector of absolute indexes.")

  .def("charge", &PLMD::pycv::PythonCVInterface::getCharge,
       "Get charge of i-th atom", py::arg("i"))

  .def("mass", &PLMD::pycv::PythonCVInterface::getMass,
       "Get mass of i-th atom", py::arg("i"))

  .def("masses",
  [](PLMD::pycv::PythonCVInterface* self) -> py::array_t<double> {
    auto nat=self->getPositions().size();
    py::array_t<double>::ShapeContainer shape({nat});
    py::array_t<double> masses(shape);
    auto retAccessor = masses.mutable_unchecked<1>();
    for (decltype(nat) i=0; i < nat; ++i) {
      retAccessor(i)
      =self->getMass(i);
    }
    return masses;
  },
  "Returns and ndarray with the masses of the atoms requested by the action")

  .def("charges",
  [](PLMD::pycv::PythonCVInterface* self) -> py::array_t<double> {
    auto nat=self->getPositions().size();
    py::array_t<double>::ShapeContainer shape({nat});
    py::array_t<double> charges(shape);
    auto retAccessor = charges.mutable_unchecked<1>();
    for (decltype(nat) i=0; i < nat; ++i) {
      retAccessor(i)
      =self->getCharge(i);
    }
    return charges;
  },
  "Returns and ndarray with the charges of the atoms requested by the action");

  /***********************************PLMD::Pbc********************************/
  py::class_<PLMD::Pbc>(m, "Pbc")
  //.def(py::init<>())
  .def("apply",[](const PLMD::Pbc* self, py::array_t<double, py::array::c_style | py::array::forcecast>& deltas) -> py::array_t<double> {
    //TODO:shape check
    //TODO: this modifies the passed deltas, so no new intialization
    //      ^this needs to be VERY clear to te user
    self->apply(PLMD::VectorView(deltas.mutable_unchecked<2>().mutable_data(0,0),deltas.shape(0)));
    return deltas;
  },
  "Apply PBC to a set of positions or distance vectors")

  //this should be optimized for speed
#define T3x3toArray( boxGetter ) \
py::array_t<double> toRet({3,3}); \
auto retAccessor = toRet.mutable_unchecked<2>(); \
PLMD::Tensor box=boxGetter; \
for(unsigned i=0;i<3;++i){ \
  for(unsigned j=0;j<3;++j) \
    retAccessor(i,j) = box(i, j); \
} \
return toRet;
  .def("getBox",[](const PLMD::Pbc* self) -> py::array_t<double> {
    T3x3toArray( self->getBox() )
  },
  "Get a numpy array of shape (3,3) with the box vectors")

  .def("getInvBox",[](const PLMD::Pbc* self) -> py::array_t<double> {
    T3x3toArray( self->getInvBox() )
  },
  "Get a numpy array of shape (3,3) with the inverted box vectors");
#undef T3x3toArray

  /******************************PLMD::NeighborList****************************/
  py::class_<PLMD::NeighborList>(m, "NeighborList")
  //.def(py::init<>())
  //https://numpy.org/doc/stable/user/basics.types.html
  //numpy.uint=unsigned long
  .def_property_readonly("size",&PLMD::NeighborList::size,"the number of pairs")

  .def("__len__", &PLMD::NeighborList::size)

  .def("getClosePairs",[](const PLMD::NeighborList* self)->py::array_t<unsigned long> {
    auto ncouples = self->size();
    py::array_t<unsigned long>::ShapeContainer shape({ncouples,2});
    py::array_t<unsigned long> couples(shape);
    auto retAccessor = couples.mutable_unchecked<2>();
    for(size_t c=0; c< ncouples; ++c) {
      retAccessor(c,0)
      = self->getClosePair(c).first;
      retAccessor(c,1) = self->getClosePair(c).second;
    }
    return couples;
  },
  "get a (NC,2) nd array with the list of couple indexes");

  /**************************PLMD::pycv::PythonFunction************************/
  py::class_<PLMD::pycv::PythonFunction>(m, "PythonFunction")
  /***************************Start of Action methods***************************/
  //usin "PLMD::pycv::PythonCVInterface::getStep" instead of the lambda gives compilation errors
  .defGetter("getStep",PLMD::pycv::PythonFunction,getStep, long int,"Returns the current step")
  .defGetter("getTime",PLMD::pycv::PythonFunction,getTime,double,"Return the present time")
  .defGetter("getTimeStep",PLMD::pycv::PythonFunction,getTimeStep,double,"Return the timestep")
  .defGetter("isRestart",PLMD::pycv::PythonFunction,getRestart,bool,"Return true if we are doing a restart")
  .defGetter("isExchangeStep",PLMD::pycv::PythonFunction,getExchangeStep,bool,"Check if we are on an exchange step")
  .def_property_readonly("label",
  [](PLMD::pycv::PythonFunction* self) -> std::string {
    return self->getLabel();
  },
  "returns the label")

  .def("log",[](PLMD::pycv::PythonFunction* self, py::object data) {
    self->log << py::str(data).cast<std::string>();
  },
  "puts a string in the PLUMED output",py::arg("s"))
  .def("lognl",[](PLMD::pycv::PythonFunction* self, py::object data) {
    self->log << py::str(data).cast<std::string>()<< "\n";
  },
  "puts a string in the PLUMED output (and appends a newline)",py::arg("s"))

  /****************************End of Action methods***************************/
  .def("argument", &PLMD::pycv::PythonFunction::
       getArgument,"Get value of the of i-th argument", py::arg("i"))

  .def("arguments", [](PLMD::pycv::PythonFunction* self) -> py::array_t<double> {
    auto nargs=self->getNumberOfArguments();
    py::array_t<double>::ShapeContainer shape({nargs});
    py::array_t<double> arguments(shape);
    for(auto i=0u; i < nargs; ++i)
      arguments.mutable_at(i) = self->getArgument(i);
    return arguments;
  }
  ,"Retuns a ndarray with the values of the arguments")

  .def_property_readonly("nargs", &PLMD::pycv::PythonFunction::
                         getNumberOfArguments,"Get the number of arguments")

  .def("difference", &PLMD::pycv::PythonFunction::
       difference,"Takes the difference taking into account pbc for argument i",
       py::arg("i"),py::arg("x"),py::arg("y"))

  .def("bringBackInPbc", &PLMD::pycv::PythonFunction::
       bringBackInPbc,"Takes one value and brings it back into the pbc of argument i",
       py::arg("i"),py::arg("x"))

  //I cannot find a way of testing properly this:
  // .def("getProjection", &PLMD::pycv::PythonFunction::
  //      getProjection,"Get the scalar product between the gradients of two variables",
  //      py::arg("i"),py::arg("j"))
  ;
}
