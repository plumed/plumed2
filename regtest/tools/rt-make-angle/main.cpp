#include "plumed/tools/Angle.h"
#include "plumed/tools/Vector.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int testCompute2Vectors90Degrees(std::ofstream &ofs) {
  ofs << "testCompute2Vectors90Degrees" << std::endl;
  Vector vector1(1.0, 0.0, 0.0);
  Vector vector2(0.0, 1.0, 0.0);
  Angle angle;
  double result = angle.compute(vector1, vector2);
  ofs << result << std::endl;
  return 0;
}

int testCompute2Vectors45Degrees(std::ofstream &ofs) {
  ofs << "testCompute2Vectors45Degrees" << std::endl;
  Vector vector1(1.0, 0.0, 0.0);
  Vector vector2(1.0, 1.0, 0.0);
  Angle angle;
  double result = angle.compute(vector1, vector2);
  ofs << result << std::endl;
  return 0;
}

int testCompute4Vectors45Degrees(std::ofstream &ofs) {
  ofs << "testCompute4Vectors45Degrees" << std::endl;
  Vector vector1(1.0, 1.0, 0.0);
  Vector vector2(1.0, 0.0, 0.0);
  Vector d1, d2;
  Angle angle;
  double result = angle.compute(vector1, vector2, d1, d2);
  ofs << result << std::endl;
  ofs << d1 << std::endl;
  ofs << d2 << std::endl;
  return 0;
}

int testCompute4VectorsParallel(std::ofstream &ofs) {
  ofs << "testCompute4VectorsParallel" << std::endl;
  Vector vector1(1.0, 1.0, 0.0);
  Vector vector2(1.0, 1.0, 0.0);
  Vector d1, d2;
  Angle angle;
  double result = angle.compute(vector1, vector2, d1, d2);
  ofs << result << std::endl;
  ofs << d1 << std::endl;
  ofs << d2 << std::endl;
  return 0;
}

int testCompute4VectorsAntiParallel(std::ofstream &ofs) {
  ofs << "testCompute4VectorsAntiParallel" << std::endl;
  Vector vector1(1.0, 1.0, 0.0);
  Vector vector2(-1.0, -1.0, 0.0);
  Vector d1, d2;
  Angle angle;
  double result = angle.compute(vector1, vector2, d1, d2);
  ofs << result << std::endl;
  ofs << d1 << std::endl;
  ofs << d2 << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testCompute2Vectors90Degrees(ofs);
  testCompute2Vectors45Degrees(ofs);
  testCompute4Vectors45Degrees(ofs);
  testCompute4VectorsParallel(ofs);
  testCompute4VectorsAntiParallel(ofs);
  return 0;
}
