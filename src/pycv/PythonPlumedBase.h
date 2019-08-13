
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
  const std::string PYTHONCV_CITATION = "(?)";

protected:
  py::module py_module;
  py::object py_fcn;

};

}
}


#endif

