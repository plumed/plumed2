#include "plumed/tools/LatticeReduction.h"
#include "plumed/tools/Tensor.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int testReduce(std::ofstream &ofs) {
  ofs << "testReduce2" << std::endl;
  Tensor tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0);
  LatticeReduction::reduce(tensor);
  ofs << tensor[0][1] << " " << tensor[0][2] << " " << tensor[0][3]
      << std::endl;
  ofs << tensor[1][1] << " " << tensor[1][2] << " " << tensor[1][3]
      << std::endl;
  ofs << tensor[2][1] << " " << tensor[2][2] << std::endl;
  return 0;
}

int testReduceSlow(std::ofstream &ofs) {
  ofs << "testReduceSlow" << std::endl;
  Tensor tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0);
  LatticeReduction::reduce(tensor);
  ofs << tensor[0][1] << " " << tensor[0][2] << " " << tensor[0][3]
      << std::endl;
  ofs << tensor[1][1] << " " << tensor[1][2] << " " << tensor[1][3]
      << std::endl;
  ofs << tensor[2][1] << " " << tensor[2][2] << std::endl;
  return 0;
}

int testReduceFast(std::ofstream &ofs) {
  ofs << "testReduceFast" << std::endl;
  Tensor tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0);
  LatticeReduction::reduce(tensor);
  ofs << tensor[0][1] << " " << tensor[0][2] << " " << tensor[0][3]
      << std::endl;
  ofs << tensor[1][1] << " " << tensor[1][2] << " " << tensor[1][3]
      << std::endl;
  ofs << tensor[2][1] << " " << tensor[2][2] << std::endl;
  return 0;
}

int testIsReduced(std::ofstream &ofs) {
  ofs << "testIsReduced" << std::endl;
  Tensor tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0);
  bool result = LatticeReduction::isReduced(tensor);
  ofs << result << std::endl;
  LatticeReduction::reduce(tensor);
  ofs << result << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testReduce(ofs);
  testReduceSlow(ofs);
  testReduceFast(ofs);
  testIsReduced(ofs);
  return 0;
}
