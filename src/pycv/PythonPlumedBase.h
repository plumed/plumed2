/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019 of Toni Giorgino

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
#ifndef __PLUMED_pycv_PythonPlumedBase_h
#define __PLUMED_pycv_PythonPlumedBase_h


#include <string>

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace PLMD {
namespace pycv {

typedef float pycv_t;		// May need to adapt to the build precision?


class PythonPlumedBase {
public:
  const std::string PYTHONCV_CITATION = "Giorgino, (2019). PYCV: a PLUMED 2 Module Enabling the Rapid Prototyping of Collective Variables in Python. Journal of Open Source Software, 4(42), 1773. doi:10.21105/joss.01773";
  const char * BIASING_DISABLED = "PYCV: Gradient was expected as a second return value but is missing. Biasing won't work\n";
  PythonPlumedBase();
  // ~PythonPlumedBase();

protected:
  py::module py_module;
  py::object py_fcn;

private:
  static int use_count;

};

}
}


#endif

