
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
  //  PythonPlumedBase();

protected:
  py::module py_module;
  py::object py_fcn;

  // We can only have one interpreter globally. This is less than ideal
  // because CVs can interfere with each other. The whole purpose of
  // this superclass is to make a singleton.
  static py::scoped_interpreter guard;

private:
  bool interpreter_initialized;

};

}
}


#endif

