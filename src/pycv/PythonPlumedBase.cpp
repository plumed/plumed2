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

  // py::scoped_interpreter PythonPlumedBase::guard{};
  
  // static py::scoped_interpreter guard{}; // start the interpreter and keep it alive

  /*
  PythonPlumedBase::PythonPlumedBase() {
    std::cout << "------ Constructor" << std::endl;
    if(!interpreter_initialized) {
      interpreter_initialized = true;
      // guard = py::scoped_interpreter{};
    }
  }
  */

}
}



