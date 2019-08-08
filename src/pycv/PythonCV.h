
#ifndef __PLUMED_pycv_PythonCV_h
#define __PLUMED_pycv_PythonCV_h


#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace PLMD {
namespace PythonCV {

  typedef float pycv_t;

  const std::string PYTHONCV_CITATION = "(you wish)";


}
}


#endif

