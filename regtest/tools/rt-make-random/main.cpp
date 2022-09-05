#include "plumed/tools/LatticeReduction.h"
#include "plumed/tools/Random.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int testConstructor(std::ofstream &ofs) {
  ofs << "testConstructor" << std::endl;
  Random random("random");
  return 0;
}

int testRandU01(std::ofstream &ofs) {
  ofs << "testRandU01" << std::endl;
  Random random("random");
  ofs << random.RandU01() << std::endl;
  random.IncreasedPrecis(true);
  ofs << random.RandU01() << std::endl;
  return 0;
}

int testU01d(std::ofstream &ofs) {
  ofs << "testU01d" << std::endl;
  Random random("random");
  ofs << random.U01d() << std::endl;
  return 0;
}

int testU01(std::ofstream &ofs) {
  ofs << "testU01" << std::endl;
  Random random("random");
  ofs << random.U01() << std::endl;
  return 0;
}

int testWriteStateFull(std::ofstream &ofs) {
  ofs << "testWriteStateFull" << std::endl;
  Random random("random");
  random.WriteStateFull(ofs);
  return 0;
}

int testStringSimple(std::ofstream &ofs) {
  ofs << "testStringSimple" << std::endl;
  Random random("random");
  std::string str = "";
  random.U01();
  random.toString(str);
  ofs << str << std::endl;
  random.fromString(str);
  return 0;
}

int testStringComplex(std::ofstream &ofs) {
  ofs << "testStringComplex" << std::endl;
  Random random("random");
  std::string str = "";
  random.toString(str);
  ofs << str << std::endl;
  random.fromString(str);
  return 0;
}

int testGaussian(std::ofstream &ofs) {
  ofs << "testGaussian" << std::endl;
  Random random("random");
  double result = random.Gaussian();
  ofs << result << std::endl;
  return 0;
}

int testShuffleSimple(std::ofstream &ofs) {
  ofs << "testShuffleSimple" << std::endl;
  Random random("random");
  std::vector<unsigned> vector = {1, 2};
  random.Shuffle(vector);
  ofs << vector[0] << " " << vector[1] << std::endl;
  return 0;
}

int testShuffleComplex(std::ofstream &ofs) {
  ofs << "testShuffleComplex" << std::endl;
  Random random("random");
  std::vector<unsigned> vector = {1, 2, 3, 4};
  random.Shuffle(vector);
  ofs << vector[0] << " " << vector[1] << " " << vector[2] << " " << vector[3] << std::endl;
  return 0;
}


int main() {
  std::ofstream ofs("output");
  testConstructor(ofs);
  testRandU01(ofs);
  testU01d(ofs);
  testU01(ofs);
  testWriteStateFull(ofs);
  testStringSimple(ofs);
  testStringComplex(ofs);
  testGaussian(ofs);
  testShuffleSimple(ofs);
  testShuffleComplex(ofs);
  return 0;
}
