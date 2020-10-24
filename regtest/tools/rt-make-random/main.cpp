#include "plumed/tools/LatticeReduction.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/Random.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int testConstructor(std::ofstream &ofs) {
  ofs << "testConstructor" << std::endl;
  new Random("random");
  return 0;
}

int testRandU01(std::ofstream &ofs) {
  ofs << "testRandU01" << std::endl;
  Random* random = new Random("random");
  ofs << random->RandU01() << std::endl;
  random->IncreasedPrecis(true);
  ofs << random->RandU01() << std::endl;
  return 0;
}

int testU01d(std::ofstream &ofs) {
  ofs << "testU01d" << std::endl;
  Random* random = new Random("random");
  ofs << random->U01d() << std::endl;
  return 0;
}

int testU01(std::ofstream &ofs) {
  ofs << "testU01" << std::endl;
  Random* random = new Random("random");
  ofs << random->U01() << std::endl;
  return 0;
}

int testWriteStateFull(std::ofstream &ofs) {
  ofs << "testWriteStateFull" << std::endl;
  Random* random = new Random("random");
  random->WriteStateFull(ofs);
  return 0;
}

int testString(std::ofstream &ofs) {
  ofs << "testString" << std::endl;
  Random* random = new Random("random");
  std::string str = "";
  random->toString(str);
  ofs << str << std::endl;
  random->fromString(str);
  return 0;
}

int testGaussian(std::ofstream &ofs) {
  ofs << "testGaussian" << std::endl;
  Random* random = new Random("random");
  double result = random->Gaussian();
  ofs << result << std::endl;
  return 0;
}

int testShuffle(std::ofstream &ofs) {
  ofs << "testShuffle" << std::endl;
  Random* random = new Random("random");
  std::vector<unsigned> vector = {1, 2};
  random->Shuffle(vector);
  ofs << vector[0] << " " << vector[1] << std::endl;
  return 0;
}

int main() {
  Stopwatch sw;
  sw.start();
  std::ofstream ofs("output");
  testConstructor(ofs);
  testRandU01(ofs);
  testU01d(ofs);
  testU01(ofs);
  testWriteStateFull(ofs);
  testString(ofs);
  testGaussian(ofs);
  testShuffle(ofs);
  sw.stop();
  std::cout << sw;
  return 0;
}
