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

namespace py=pybind11;


PYBIND11_EMBEDDED_MODULE(plumedCommunications, m) {
  py::class_<PLMD::pycv::PythonCVInterface>(m, "PythonCVInterface")
  .def("getStep", [](PLMD::pycv::PythonCVInterface* self) {
    return self->getStep();
  },"Returns the current step")
  .def("getPosition",&PLMD::pycv::PythonCVInterface:: getPosition,
       "Returns the vector with the position of the \"i\"th"
       " atom requested by the action",py::arg("i"))
  .def("getPositionArray",
      [](PLMD::pycv::PythonCVInterface* self, int i) -> py::array_t<double>{
        py::array_t<double>::ShapeContainer shape({3});
        py::array_t<double> atom(shape,&self->getPosition(i)[0]);
        return atom;
      },
       "Returns an ndarray with the position of the \"i\"th"
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
  .def("getPbc",&PLMD::pycv::PythonCVInterface::getPbc)
  ;
  py::class_<PLMD::Pbc>(m, "PLMDPbc")
  //.def(py::init<>())
  .def("apply",[](PLMD::Pbc* self, py::array_t<double>& deltas) -> py::array_t<double>{
    //TODO:shape check
    auto accessor = deltas.unchecked<2>(); 
    auto nat=deltas.shape(0);
    py::array_t<double> toRet(py::array_t<double>::ShapeContainer({nat,3ul}));
    auto retAccessor = toRet.mutable_unchecked<2>(); 
    for (auto i=0u; i < nat;++i){
      auto t= self->distance(PLMD::Vector(0.0,0.0,0.0),
        //I think this may be slow, but thys is a prototype
        PLMD::Vector(accessor(i,0),accessor(i,1),accessor(i,2)),nullptr);
        retAccessor(i,0) = t[0];
        retAccessor(i,1) = t[1];
        retAccessor(i,2) = t[2];
    }
    return toRet;
  },
  "Apply PBC to a set of positions or distance vectors")
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
