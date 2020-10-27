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

int testGetDx(std::ofstream &ofs) {
  ofs << "testGetDx" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<double> property = gridBase->getDx();
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ';
  }
  ofs << std::endl;
  return 0;
}

int testGetBinVolume(std::ofstream &ofs) {
  ofs << "testGetBinVolume" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  double binVolume = gridBase->getBinVolume();
  ofs << binVolume << std::endl;
  return 0;
}

int testGetNbin(std::ofstream &ofs) {
  ofs << "testGetNbin" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<unsigned> property = gridBase->getNbin();
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ';
  }
  ofs << std::endl;
  return 0;
}

int testGetIsPeriodic(std::ofstream &ofs) {
  ofs << "testGetIsPeriodic" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<bool> property = gridBase->getIsPeriodic();
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ';
  }
  ofs << std::endl;
  return 0;
}

int testGetDimension(std::ofstream &ofs) {
  ofs << "testGetDimension" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  unsigned dimension = gridBase->getDimension();
  ofs << dimension << std::endl;
  return 0;
}

int testGetArgNames(std::ofstream &ofs) {
  ofs << "testGetArgNames" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<std::string> property = gridBase->getArgNames();
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ';
  }
  ofs << std::endl;
  return 0;
}

int testHasDerivatives(std::ofstream &ofs) {
  ofs << "testHasDerivatives" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  unsigned property = gridBase->hasDerivatives();
  ofs << property << std::endl;
  return 0;
}

int testGetIndiciesByIndex(std::ofstream &ofs) {
  ofs << "testGetIndiciesByIndex" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<unsigned> property = gridBase->getIndices(4);
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ';
  }
  ofs << std::endl;
  return 0;
}

int testGetIndiciesByPointVector(std::ofstream &ofs) {
  ofs << "testGetIndiciesByPointVector" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  const std::vector<double> vector = {1.32};
  std::vector<unsigned> property = gridBase->getIndices(vector);
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
  testGetDx(ofs);
  testGetBinVolume(ofs);
  testGetNbin(ofs);
  testGetIsPeriodic(ofs);
  testGetDimension(ofs);
  testGetArgNames(ofs);
  testHasDerivatives(ofs);
  testGetIndiciesByIndex(ofs);
  testGetIndiciesByPointVector(ofs);
  return 0;
}
