#include "plumed/tools/LatticeReduction.h"
#include "plumed/tools/Tensor.h"
#include <fstream>
#include <iostream>
#include <memory>

using namespace PLMD;

int testReduce(std::ofstream &ofs) {
  ofs << "testReduce2" << std::endl;
  std::unique_ptr<Tensor> tensor (new Tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0));
  LatticeReduction::reduce(*tensor);
  ofs << tensor->getRow(0)[0] << " " << tensor->getRow(0)[1] << " " << tensor->getRow(0)[2] << std::endl;
  ofs << tensor->getRow(1)[0] << " " << tensor->getRow(1)[1] << " " << tensor->getRow(1)[2] << std::endl;
  ofs << tensor->getRow(2)[0] << " " << tensor->getRow(2)[1] << " " << tensor->getRow(2)[2] << std::endl;
  return 0;
}

int testReduceSlow(std::ofstream &ofs) {
  ofs << "testReduceSlow" << std::endl;
  std::unique_ptr<Tensor> tensor (new Tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0));
  LatticeReduction::reduceSlow(*tensor);
  ofs << tensor->getRow(0)[0] << " " << tensor->getRow(0)[1] << " " << tensor->getRow(0)[2] << std::endl;
  ofs << tensor->getRow(1)[0] << " " << tensor->getRow(1)[1] << " " << tensor->getRow(1)[2] << std::endl;
  ofs << tensor->getRow(2)[0] << " " << tensor->getRow(2)[1] << " " << tensor->getRow(2)[2] << std::endl;
  return 0;
}

int testReduceFast(std::ofstream &ofs) {
  ofs << "testReduceFast" << std::endl;
  std::unique_ptr<Tensor> tensor (new Tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0));
  LatticeReduction::reduceFast(*tensor);
  ofs << tensor->getRow(0)[0] << " " << tensor->getRow(0)[1] << " " << tensor->getRow(0)[2] << std::endl;
  ofs << tensor->getRow(1)[0] << " " << tensor->getRow(1)[1] << " " << tensor->getRow(1)[2] << std::endl;
  ofs << tensor->getRow(2)[0] << " " << tensor->getRow(2)[1] << " " << tensor->getRow(2)[2] << std::endl;
  return 0;
}

int testIsReduced(std::ofstream &ofs) {
  ofs << "testIsReduced" << std::endl;
  std::unique_ptr<Tensor> tensor (new Tensor(1.0, 2.0, 3.0, 5.0, 4.0, 3.0, 10.0, 8.0, 2.0));
  bool result1 = LatticeReduction::isReduced(*tensor);
  ofs << result1 << std::endl;
  LatticeReduction::reduce(*tensor);
  bool result2 = LatticeReduction::isReduced(*tensor);
  ofs << result2 << std::endl;
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
