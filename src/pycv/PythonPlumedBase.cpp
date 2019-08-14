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



namespace py = pybind11;


namespace PLMD {
namespace pycv {

// Unfortunately we can only have one interpreter globally. This is
// less than ideal because CVs can interfere with each other. The
// whole purpose of this superclass is to make a singleton with it.
static py::scoped_interpreter guard{}; // start the interpreter and keep it alive

}
}



