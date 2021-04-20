#include "plumed/tools/h36.h"
#include <fstream>
#include <iostream>

using namespace PLMD::h36;

int testHy36encode(std::ofstream &ofs) {
  ofs << "testHy36encode" << std::endl;
  char result[4 + 1];
  hy36encode(4, 12345, result);
  ofs << "result= " << result << std::endl;
  return 0;
}

int testHy36decodeWidth4(std::ofstream &ofs) {
  ofs << "testHy36decodeWidth4" << std::endl;
  int result;
  hy36decode(4, "A1T5", 4, &result);
  ofs << "result= " << result << std::endl;
  return 0;
}

int testHy36decodeWidth5(std::ofstream &ofs) {
  ofs << "testHy36decodeWidth5" << std::endl;
  int result;
  hy36decode(5, "A1T5", 5, &result);
  ofs << "result= " << result << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testHy36encode(ofs);
  testHy36decodeWidth4(ofs);
  testHy36decodeWidth5(ofs);
  return 0;
}
