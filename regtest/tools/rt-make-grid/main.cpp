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
  delete value;
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
    ofs << i << " " << property[i] << ' ' << std::endl;
  }
  return 0;
}

int testGetIndiciesByPointVector(std::ofstream &ofs) {
  ofs << "testGetIndiciesByPointVector" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  const std::vector<double> vector = {1.32};
  std::vector<unsigned> property = gridBase->getIndices(vector);
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ' << std::endl;
  }
  return 0;
}

int testGetPointByIndex(std::ofstream &ofs) {
  ofs << "testGetPointByIndex" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<double> property = gridBase->getPoint(0);
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ' << std::endl;
  }
  return 0;
}

int testGetNearestNeighbors(std::ofstream &ofs) {
  ofs << "testGetNearestNeighbors" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  std::vector<GridBase::index_t> property = gridBase->getNearestNeighbors(1);
  for (int i = 0; i < property.size(); ++i) {
    ofs << i << " " << property[i] << ' ' << std::endl;
  }
  return 0;
}

int testGetSize(std::ofstream &ofs) {
  ofs << "testGetSize" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  unsigned property = gridBase->getSize();
  ofs << property << std::endl;
  return 0;
}

int testGetValue(std::ofstream &ofs) {
  ofs << "testGetValue" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  double property = gridBase->getValue(0);
  ofs << property << std::endl;
  return 0;
}

int testAddValue(std::ofstream &ofs) {
  ofs << "testAddValue" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  ofs << gridBase->getValue(0) << std::endl;
  gridBase->addValue(0, 1.0);
  ofs << gridBase->getValue(0) << std::endl;
  return 0;
}

int testOutputFormat(std::ofstream &ofs) {
  ofs << "testOutputFormat" << std::endl;
  std::unique_ptr<GridBase> gridBase = createGrid();
  gridBase->setOutputFmt("%14.9f");
  gridBase->resetToDefaultOutputFmt();
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
  testGetPointByIndex(ofs);
  testGetNearestNeighbors(ofs);
  testGetSize(ofs);
  testGetValue(ofs);
  testAddValue(ofs);
  testOutputFormat(ofs);
  return 0;
}
