#include "plumed/core/Action.h"
#include "plumed/core/ActionWithValue.h"
#include "plumed/core/PlumedMain.h"
#include "plumed/core/Value.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace PLMD;

int testConstructor(std::ofstream &ofs) {
  ofs << "testConstructor" << std::endl;
  Value *value = new Value();
  ofs << value->getMaxMinusMin() << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testConstructor(ofs);
  return 0;
}
