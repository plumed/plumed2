#include <iostream>
#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace py::literals;

// https://pybind11.readthedocs.io/en/stable/advanced/embedding.html
// https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html

py::scoped_interpreter guard{}; // start the interpreter and keep it alive

void main1() {
    py::print("Hello, World!"); // use the Python API

    auto locals = py::dict("name"_a="World", "number"_a=42);
    py::exec(R"(
        message = "Hello, {name}! The answer is {number}".format(**locals())
    )", py::globals(), locals);

    auto message = locals["message"].cast<std::string>();
    std::cout << message;
}


void main2() {
  const int N=3, M=5;

  py::array_t<float, py::array::c_style> arr({ N, M });
  float * ptr= (float*) arr.request().ptr;
  for(int i=0; i<N; i++) {
    for(int j=0; j<M; j++) {
      // arr.data(i,j) = i*100+j;
      ptr[i*M+j]  = i*100+j;
    }
  }

  py::module mod = py::module::import("mybind");
  auto result = mod.attr("callme")(arr);

  for(int i=0; i<N; i++) {
    for(int j=0; j<M; j++) {
      // arr.data(i,j) = i*100+j;
      ptr[i*M+j]  = i*200+j;
    }
  }

  result = mod.attr("callme")(arr);
  py::array_t<float> result_a(result);

  std::cout << ptr[0] << std::endl;

  auto r= result_a.unchecked<2>();
  std::cout << r(0,0) << std::endl;

  
}



int main() {
  main2();
}
