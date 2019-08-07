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

  // py::array_t<float, py::array::c_style> arr({ N, M });
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
      arr.mutable_at(i,j) = i*200+j;
      // ptr[i*M+j]  = i*200+j;
    }
  }


  result = mod.attr("callme")(arr);
  py::array_t<float> result_a(result);

  std::cout << "T1 " << ptr[0] << std::endl;

  auto r= result_a.unchecked<2>();
  std::cout << "T2 " << r(0,0) << std::endl;

  std::cout << "T3 " << result_a.at(0,0) << std::endl;

  std::cout << "T4 " << arr.size() << arr.ndim() << arr.shape(0) << arr.shape(1) << std::endl;
 

  
}


void main3() {
    py::module mod = py::module::import("mybind");
    py::list r = mod.attr("main3")();
    // py::list rl=r.cast<py::list>();
    int r1=r[1].cast<int>();
    std::cout << r1 << std::endl;
}


int main() {
  main2();
}
