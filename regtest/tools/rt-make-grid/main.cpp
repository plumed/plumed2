#include "plumed/core/Value.h"
#include "plumed/tools/File.h"
#include "plumed/tools/Grid.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int testCreate(std::ofstream &ofs) {
  ofs << "testCreate" << std::endl;
  IFile gridfile;
  gridfile.open("./input-grid.data");

  Value *value = new Value("d1");
  value->setNotPeriodic();
  const std::vector<Value *> grid_args_ = {value};
  std::unique_ptr<GridBase> gridBase = GridBase::create(
      "external.bias", grid_args_, gridfile, false, false, false);
  ofs << "dimensions= " << gridBase->getDimension() << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testCreate(ofs);
  return 0;
}
