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

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

#include <iostream>


namespace py = pybind11;


namespace PLMD {
namespace pycv {

int PythonPlumedBase::use_count=0;

// We can only have one interpreter globally. This is less than ideal
// because CVs can interfere with each other. The whole purpose of
// this superclass is to make a singleton.
// https://pybind11.readthedocs.io/en/stable/reference.html#_CPPv422initialize_interpreterb

PythonPlumedBase::PythonPlumedBase() {
  if(use_count==0) {
    // std::cout << "------ init" << std::endl;
    py::initialize_interpreter();
  } else {
    // std::cout << "------ reusing" << std::endl;
  }
  use_count++;
}

// Finalization is tricky, because it should happen AFTER the
// destruction of the derived classes (which contain py::
// objects). Not doing it.


}
}



