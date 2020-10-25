#include "plumed/core/Value.h"
#include "plumed/tools/File.h"
#include "plumed/tools/Grid.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

std::unique_ptr<GridBase> createGrid() {
  IFile gridfile;
  gridfile.open("./input-grid.data");

  Value *value = new Value("d1");
  value->setNotPeriodic();
  const std::vector<Value *> grid_args_ = {value};
  std::unique_ptr<GridBase> gridBase = GridBase::create(
      "external.bias", grid_args_, gridfile, false, false, false); 
  return gridBase;
}

int testCreate(std::ofstream &ofs) {
  ofs << "testCreate" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  ofs << "dimensions= " << gridBase->getDimension() << std::endl;
  return 0;
}

int testGetMax(std::ofstream &ofs) {
  ofs << "testGetMax" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<std::string> property = gridBase->getMax();
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ';
  }
  ofs << std::endl;
  return 0;
}

int testGetMin(std::ofstream &ofs) {
  ofs << "testGetMin" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<std::string> property = gridBase->getMin();
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ';
  }
  ofs << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testCreate(ofs);
  testGetMax(ofs);
  testGetMin(ofs);
  return 0;
}
