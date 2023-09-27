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
#include <pybind11/operators.h>

#include "tools/Vector.h"

#include "PythonCVInterface.h"
#include "PythonPlumedBase.h"

namespace py=pybind11;
using PLMD::pycv::pycv_t;

PYBIND11_EMBEDDED_MODULE(plumedCommunications, m) {
  py::class_<PLMD::pycv::PythonCVInterface>(m, "PythonCVInterface")
  //usin "PLMD::pycv::PythonCVInterface::getStep" instead of the lambda gives compilation errors
  .def("getStep",
       [](PLMD::pycv::PythonCVInterface* self) -> long int{return self->getStep();},
       "Returns the current step")
  .def("getPosition",
  [](PLMD::pycv::PythonCVInterface* self, int i) -> py::array_t<double> {
    py::array_t<double>::ShapeContainer shape({3});
    py::array_t<double> atom(shape,&self->getPosition(i)[0]);
    return atom;
  },
  "Returns an ndarray with the position of the \"i\"th"
  " atom requested by the action",py::arg("i"))
  .def("getPositions",[](PLMD::pycv::PythonCVInterface* self) -> py::array_t<double> {
    auto nat=self->getPositions().size();
    py::array_t<double>::ShapeContainer shape({nat,3});
    py::array_t<double> atomList(shape,
                                 &self->getPositions()[0][0]
                                );
    return atomList;
  },"Returns a numpy.array that contaisn the atomic positions of the atoms requested by the action")
  .def("getPbc",&PLMD::pycv::PythonCVInterface::getPbc,
       "returns an interface to the current pbcs")
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
                             PLMD::Vector(accessor(i,0),accessor(i,1),accessor(i,2)),nullptr);
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
}
