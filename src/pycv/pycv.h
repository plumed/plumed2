
#ifndef __PLUMED_pycv_pycv_h
#define __PLUMED_pycv_pycv_h


#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace PLMD {
namespace pycv {

  const std::string PYTHONCV_CITATION = "(?)";

  
  // Unfortunately we can only have one interpreter globally. This is
  // less than ideal because CVs can interfere with each other.
  static py::scoped_interpreter guard{}; // start the interpreter and keep it alive 

}
}


#endif

