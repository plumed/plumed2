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

#include "tools/Vector.h"
#include "tools/NeighborList.h"

#include "PythonCVInterface.h"

namespace py=pybind11;
using PLMD::pycv::pycv_t;

PYBIND11_MAKE_OPAQUE(std::vector<PLMD::AtomNumber>)

#define defGetter(pyfun, cppfun, type, description) \
  def(pyfun, [](PLMD::pycv::PythonCVInterface* self) -> type{ \
    return self->cppfun(); }, description)

PYBIND11_EMBEDDED_MODULE(plumedCommunications, m) {
  py::module_ defaults = m.def_submodule("defaults", "Submodule with the default definitions");
  defaults.attr("COMPONENT") = py::dict(py::arg("period")=py::none(),py::arg("derivative")=true);
  defaults.attr("COMPONENT_NODEV") = py::dict(py::arg("period")=py::none(),py::arg("derivative")=false);
  py::bind_vector<std::vector<PLMD::AtomNumber>>(m, "VectorAtomNumber");
  py::class_<PLMD::pycv::PythonCVInterface>(m, "PythonCVInterface")
  .def_readwrite("data",&PLMD::pycv::PythonCVInterface::dataContainer,"Return an accessible dictionary that persist along all the simulation")
  //usin "PLMD::pycv::PythonCVInterface::getStep" instead of the lambda gives compilation errors
  .defGetter("getStep",getStep, long int,"Returns the current step")
  .defGetter("getTime",getTime,double,"Return the present time")
  .defGetter("getTimeStep",getTimeStep,double,"Return the timestep")
  .defGetter("isRestart",getRestart,bool,"Return true if we are doing a restart")
  .defGetter("isExchangeStep",getExchangeStep,bool,"Check if we are on an exchange step")
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
  .def_property_readonly("label",
  [](PLMD::pycv::PythonCVInterface* self) -> std::string {
    return self->getLabel();
  },
  "returns the label")
  //here I can use &PLMD::pycv::PythonCVInterface::makeWhole because is not in Action
  // that is inherithed both by withValue and Atomistic
  .def("makeWhole",&PLMD::pycv::PythonCVInterface::makeWhole,"Make atoms whole, assuming they are in the proper order")
  .def("getPositions",[](PLMD::pycv::PythonCVInterface* self) -> py::array_t<double> {
    auto nat=self->getPositions().size();
    py::array_t<double>::ShapeContainer shape({nat,3});
    py::array_t<double> atomList(shape,
                                 &self->getPositions()[0][0]
                                );
    return atomList;
  },
  "Returns a numpy.array that contains the atomic positions of the atoms requested by the action")
  .def("log",[](PLMD::pycv::PythonCVInterface* self, py::object data) {
    self->log << py::str(data).cast<std::string>();
  },
  "put a string in the PLUMED output",py::arg("s"))
  .def("lognl",[](PLMD::pycv::PythonCVInterface* self, py::object data) {
    self->log << py::str(data).cast<std::string>()<< "\n";
  },
  "put a string in the PLUMED output (it appends a newline)",py::arg("s"))
  .def("getPbc",&PLMD::pycv::PythonCVInterface::getPbc,
       "returns an interface to the current pbcs")
  .def("getNeighbourList",&PLMD::pycv::PythonCVInterface::getNL,
       "returns an interface to the current Neighborlist")
//https://pybind11.readthedocs.io/en/stable/advanced/functions.html#return-value-policies
  /*.def("getAbsoluteIndexes", &PLMD::pycv::PythonCVInterface::getAbsoluteIndexes ,
  "Get the vector of absolute indexes.",
  py::return_value_policy::reference)*/
  .def_property_readonly("absoluteIndexes",
                         &PLMD::pycv::PythonCVInterface::getAbsoluteIndexes,
                         "Get the vector of absolute indexes.",
                         py::return_value_policy::reference_internal
                        )
  .def("charge", &PLMD::pycv::PythonCVInterface::getCharge,
       "Get charge of i-th atom", py::arg("i")
      )
  .def("mass", &PLMD::pycv::PythonCVInterface::getMass,
       "Get mass of i-th atom", py::arg("i")
      )
  .def("masses",
  [](PLMD::pycv::PythonCVInterface* self) -> py::array_t<double> {
    auto nat=self->getPositions().size();
    py::array_t<double>::ShapeContainer shape({nat});
    py::array_t<double> masses(shape);
    auto retAccessor = masses.mutable_unchecked<1>();
    for (auto i=0u; i < nat; ++i) {
      retAccessor(i)=self->getMass(i);
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
    for (auto i=0u; i < nat; ++i) {
      retAccessor(i)=self->getCharge(i);
    }
    return charges;
  },
  "Returns and ndarray with the charges of the atoms requested by the action")
  ;
  py::class_<PLMD::Pbc>(m, "PLMDPbc")
  //.def(py::init<>())
  .def("apply",[](const PLMD::Pbc* self, py::array_t<double>& deltas) -> py::array_t<double> {
    //TODO:shape check
    //TODO: this may be set up to modify the passed deltas, so no new intialization
    //      ^this needs to be VERY clear to te user
    auto accessor = deltas.unchecked<2>();
    auto nat=deltas.shape(0);
    py::array_t<double> toRet({nat,deltas.shape(1)});
    auto retAccessor = toRet.mutable_unchecked<2>();
    for (auto i=0u; i < nat; ++i) {
      auto t= self->distance(PLMD::Vector(0.0,0.0,0.0),
                             //I think this may be slow, but serves as a demonstration as a base for the tests
                             PLMD::Vector(accessor(i,0),accessor(i,1),accessor(i,2)),
                             nullptr);
      retAccessor(i,0) = t[0];
      retAccessor(i,1) = t[1];
      retAccessor(i,2) = t[2];
    }
    return toRet;
  },
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
  "Apply PBC to a set of positions or distance vectors")
  .def("getBox",[](const PLMD::Pbc* self) -> py::array_t<double> {
    T3x3toArray( self->getBox() )
    /*
    py::array_t<double> toRet({3,3});
    auto retAccessor = toRet.mutable_unchecked<2>();
    PLMD::Tensor box=self->getBox();
    //TODO: optimize for speed:
    for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j)
        retAccessor(i,j) = box(i,j);
    }
    return toRet;*/
  },
  "Get a numpy array of shape (3,3) with the box vectors")
  .def("getInvBox",[](const PLMD::Pbc* self) -> py::array_t<double> {
    T3x3toArray( self->getInvBox() )
    /*
    py::array_t<double> toRet({3,3});
    auto retAccessor = toRet.mutable_unchecked<2>();
    PLMD::Tensor box=self->getInvBox();
    //TODO: optimize for speed:
    for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j)
        retAccessor(i,j) = box(i,j);
    }
    return toRet;
    */
  },
  "Get a numpy array of shape (3,3) with the inverted box vectors")
#undef T3x3toArray
  ;
  py::class_<PLMD::NeighborList>(m, "NeighborList")
  //.def(py::init<>())
  //https://numpy.org/doc/stable/user/basics.types.html
  //numpy.uint=unsigned long
  .def("size",&PLMD::NeighborList::size,"return the number of pairs")
  .def("__len__", &PLMD::NeighborList::size)
  .def("getClosePairs",[](const PLMD::NeighborList* self)->py::array_t<unsigned long> {
    auto ncouples = self->size();
    py::array_t<unsigned long>::ShapeContainer shape({ncouples,2});
    py::array_t<unsigned long> couples(shape);
    auto retAccessor = couples.mutable_unchecked<2>();
    for(size_t c=0; c< ncouples; ++c) {
      retAccessor(c,0) = self->getClosePair(c).first;
      retAccessor(c,1) = self->getClosePair(c).second;
    }
    return couples;
  },

  "get a (NC,2) nd array with the list of couple indexes")
  ;
  py::class_<PLMD::AtomNumber>(m, "AtomNumber")
  .def(py::init<>())
//
  .def_property_readonly("index",
                         static_cast<unsigned(PLMD::AtomNumber::*)()const>(&PLMD::AtomNumber::index),
                         "The index number.")
  .def_property_readonly("serial",
                         static_cast<unsigned(PLMD::AtomNumber::*)()const>(&PLMD::AtomNumber::serial),
                         "The index number.")
  ;

}
