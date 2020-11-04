#include "plumed/core/Action.h"
#include "plumed/core/ActionWithValue.h"
#include "plumed/core/PlumedMain.h"
#include "plumed/core/Value.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace PLMD;

int testConstructor(std::ofstream &ofs) {
  ofs << "testConstructor" << std::endl;
  std::unique_ptr<Value> value (new Value());
  ofs << value->get() << std::endl;
  return 0;
}

int testConstructorWithName(std::ofstream &ofs) {
  ofs << "testConstructorWithName" << std::endl;
  std::unique_ptr<Value> value (new Value());
  ofs << value->get() << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testConstructor(ofs);
  testConstructorWithName(ofs);
  return 0;
}
